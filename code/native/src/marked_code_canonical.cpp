#include "utsp/marked_code_canonical.hpp"

#include "utsp/tensor_canonical.hpp"
#include "utsp/tensor_primitives.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <iterator>
#include <optional>
#include <set>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef UTSP_HAVE_BLISS
#include <bliss/graph.hh>
#include <bliss/stats.hh>
#endif

namespace utsp {
namespace {

[[nodiscard]] std::uint32_t pivot_of(Mask value) {
  if (value == 0) throw std::invalid_argument("zero vector has no pivot");
  return 63U - static_cast<std::uint32_t>(std::countl_zero(value));
}

[[nodiscard]] std::optional<Mask> coordinates_in_basis(
    Mask vector, std::span<const Mask> basis) {
  if (basis.size() > 64) {
    throw std::invalid_argument("coordinate basis exceeds 64 vectors");
  }
  std::array<Mask, 64> pivot_rows{};
  std::array<Mask, 64> pivot_coordinates{};
  for (std::size_t index = 0; index < basis.size(); ++index) {
    Mask row = basis[index];
    Mask coordinates = Mask{1} << index;
    while (row != 0) {
      const auto pivot = pivot_of(row);
      if (pivot_rows[pivot] != 0) {
        row ^= pivot_rows[pivot];
        coordinates ^= pivot_coordinates[pivot];
        continue;
      }
      pivot_rows[pivot] = row;
      pivot_coordinates[pivot] = coordinates;
      break;
    }
    if (row == 0) {
      throw std::invalid_argument("coordinate basis is linearly dependent");
    }
  }
  Mask coordinates = 0;
  while (vector != 0) {
    const auto pivot = pivot_of(vector);
    if (pivot_rows[pivot] == 0) return std::nullopt;
    vector ^= pivot_rows[pivot];
    coordinates ^= pivot_coordinates[pivot];
  }
  return coordinates;
}

[[nodiscard]] std::vector<Mask> physical_label_rows(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows) {
  const auto canonical_rows = rref_basis(quotient_rows);
  std::vector<Mask> result;
  result.reserve(canonical_rows.size());
  for (const auto row : canonical_rows) {
    if ((row >> label_space.quotient_dimension()) != 0) {
      throw std::invalid_argument("quotient row lies outside label space");
    }
    result.push_back(linear_combination(row, label_space.quotient_basis));
  }
  return result;
}

[[nodiscard]] std::uint32_t cubic_pair_bit(
    std::uint32_t q, std::uint32_t left, std::uint32_t right) {
  if (left >= right || right >= q) {
    throw std::invalid_argument("invalid cubic pair coefficient");
  }
  std::uint32_t bit = q;
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second, ++bit) {
      if (first == left && second == right) return bit;
    }
  }
  throw std::logic_error("cubic pair coefficient was not found");
}

[[nodiscard]] std::uint32_t cubic_triple_bit(
    std::uint32_t q,
    std::uint32_t first,
    std::uint32_t second,
    std::uint32_t third) {
  if (first >= second || second >= third || third >= q) {
    throw std::invalid_argument("invalid cubic triple coefficient");
  }
  std::uint32_t bit = q + q * (q - 1) / 2;
  for (std::uint32_t left = 0; left < q; ++left) {
    for (std::uint32_t middle = left + 1; middle < q; ++middle) {
      for (std::uint32_t right = middle + 1; right < q; ++right, ++bit) {
        if (left == first && middle == second && right == third) return bit;
      }
    }
  }
  throw std::logic_error("cubic triple coefficient was not found");
}

class CubicExtensionSignature {
 public:
  explicit CubicExtensionSignature(std::span<const Mask> parent_rows)
      : parent_rows_(parent_rows.begin(), parent_rows.end()),
        child_dimension_(
            static_cast<std::uint32_t>(parent_rows.size() + 1)) {
    if (child_dimension_ > 7) return;
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      if (std::popcount(parent_rows_[first]) & 1U) {
        embedded_parent_signature_ |= std::uint64_t{1} << first;
      }
    }
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      for (std::uint32_t second = first + 1;
           second < parent_rows_.size(); ++second) {
        if (std::popcount(parent_rows_[first] & parent_rows_[second]) & 1U) {
          embedded_parent_signature_ |=
              std::uint64_t{1}
              << cubic_pair_bit(child_dimension_, first, second);
        }
      }
    }
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      for (std::uint32_t second = first + 1;
           second < parent_rows_.size(); ++second) {
        for (std::uint32_t third = second + 1;
             third < parent_rows_.size(); ++third) {
          if (std::popcount(
                  parent_rows_[first] & parent_rows_[second] &
                  parent_rows_[third]) &
              1U) {
            embedded_parent_signature_ |=
                std::uint64_t{1}
                << cubic_triple_bit(
                       child_dimension_, first, second, third);
          }
        }
      }
    }
    const auto extension = child_dimension_ - 1;
    extension_pair_bits_.reserve(parent_rows_.size());
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      extension_pair_bits_.push_back(
          cubic_pair_bit(child_dimension_, first, extension));
    }
    extension_triples_.reserve(
        parent_rows_.size() * (parent_rows_.size() - 1) / 2);
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      for (std::uint32_t second = first + 1;
           second < parent_rows_.size(); ++second) {
        extension_triples_.push_back(ExtensionTriple{
            first,
            second,
            cubic_triple_bit(
                child_dimension_, first, second, extension),
        });
      }
    }
  }

  [[nodiscard]] bool available() const { return child_dimension_ <= 7; }

  [[nodiscard]] std::uint64_t extend(Mask extension_row) const {
    if (!available()) {
      throw std::logic_error("packed cubic extension exceeds q=7");
    }
    auto signature = embedded_parent_signature_;
    const auto extension = child_dimension_ - 1;
    if (std::popcount(extension_row) & 1U) {
      signature |= std::uint64_t{1} << extension;
    }
    for (std::uint32_t first = 0; first < parent_rows_.size(); ++first) {
      if (std::popcount(parent_rows_[first] & extension_row) & 1U) {
        signature |= std::uint64_t{1}
                     << extension_pair_bits_[first];
      }
    }
    for (const auto& triple : extension_triples_) {
      if (std::popcount(
              parent_rows_[triple.first] & parent_rows_[triple.second] &
              extension_row) &
          1U) {
        signature |= std::uint64_t{1} << triple.bit;
      }
    }
    return signature;
  }

 private:
  struct ExtensionTriple {
    std::uint32_t first;
    std::uint32_t second;
    std::uint32_t bit;
  };

  std::vector<Mask> parent_rows_;
  std::vector<std::uint32_t> extension_pair_bits_;
  std::vector<ExtensionTriple> extension_triples_;
  std::uint32_t child_dimension_ = 0;
  std::uint64_t embedded_parent_signature_ = 0;
};

#ifdef UTSP_HAVE_BLISS

[[nodiscard]] Mask permute_row(
    Mask row, std::span<const std::uint32_t> coordinate_images) {
  Mask image = 0;
  while (row != 0) {
    const auto coordinate = static_cast<std::size_t>(std::countr_zero(row));
    if (coordinate >= coordinate_images.size()) {
      throw std::invalid_argument("row lies outside coordinate permutation");
    }
    image |= Mask{1} << coordinate_images[coordinate];
    row &= row - 1;
  }
  return image;
}

[[nodiscard]] std::vector<Mask> induce_quotient_map(
    const LabelSpace& label_space,
    std::span<const std::uint32_t> coordinate_images) {
  if (coordinate_images.size() != label_space.points.size()) {
    throw std::invalid_argument("coordinate automorphism has the wrong degree");
  }
  std::vector<Mask> full_basis = label_space.stabilizer_basis;
  full_basis.insert(
      full_basis.end(),
      label_space.quotient_basis.begin(),
      label_space.quotient_basis.end());
  const auto stabilizer_dimension = label_space.stabilizer_basis.size();
  for (const auto row : label_space.stabilizer_basis) {
    const auto coordinates = coordinates_in_basis(
        permute_row(row, coordinate_images), full_basis);
    if (!coordinates || (*coordinates >> stabilizer_dimension) != 0) {
      throw std::logic_error(
          "marked-code automorphism did not preserve the stabilizer subcode");
    }
  }
  std::vector<Mask> quotient_images;
  quotient_images.reserve(label_space.quotient_dimension());
  for (const auto row : label_space.quotient_basis) {
    const auto coordinates = coordinates_in_basis(
        permute_row(row, coordinate_images), full_basis);
    if (!coordinates) {
      throw std::logic_error(
          "marked-code automorphism did not preserve the label space");
    }
    quotient_images.push_back(*coordinates >> stabilizer_dimension);
  }
  if (gf2_rank(quotient_images) != label_space.quotient_dimension()) {
    throw std::logic_error("induced quotient automorphism is singular");
  }
  return quotient_images;
}

#endif

[[nodiscard]] std::vector<Mask> all_codewords(
    std::span<const Mask> basis) {
  if (basis.size() >= 24) {
    throw std::invalid_argument(
        "marked-code graph is limited to code dimension below 24");
  }
  const auto count = std::size_t{1} << basis.size();
  std::vector<Mask> result;
  result.reserve(count);
  for (std::size_t coefficients = 0; coefficients < count; ++coefficients) {
    result.push_back(linear_combination(coefficients, basis));
  }
  return result;
}

[[nodiscard]] std::vector<Mask> marked_code_invariant(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows) {
  const auto logical_rows = physical_label_rows(label_space, quotient_rows);
  const auto stabilizer_words = all_codewords(label_space.stabilizer_basis);
  if (logical_rows.size() >= 20) {
    throw std::invalid_argument("marked-code invariant requires q < 20");
  }
  std::vector<std::vector<Mask>> coset_profiles;
  const auto logical_count = std::size_t{1} << logical_rows.size();
  coset_profiles.reserve(logical_count - 1);
  for (std::size_t coefficients = 1; coefficients < logical_count;
       ++coefficients) {
    const auto representative = linear_combination(coefficients, logical_rows);
    std::vector<Mask> profile(label_space.points.size() + 1, 0);
    for (const auto stabilizer : stabilizer_words) {
      ++profile[std::popcount(representative ^ stabilizer)];
    }
    // At the first weight where two histograms differ, the expanded sorted
    // profile with more entries at that weight is lexicographically smaller.
    // Complementing each count preserves the exact old ordering compactly.
    for (auto& count : profile) {
      count = stabilizer_words.size() - count;
    }
    coset_profiles.push_back(std::move(profile));
  }
  std::sort(coset_profiles.begin(), coset_profiles.end());
  std::vector<Mask> result;
  result.reserve(
      1 + (logical_count - 1) * (label_space.points.size() + 1));
  result.push_back(
      static_cast<Mask>(label_space.points.size()) |
      (static_cast<Mask>(label_space.stabilizer_basis.size()) << 8U) |
      (static_cast<Mask>(logical_rows.size()) << 16U));
  for (const auto& profile : coset_profiles) {
    result.insert(result.end(), profile.begin(), profile.end());
  }
  return result;
}

class HyperplaneInvariantCache {
 public:
  HyperplaneInvariantCache(
      const LabelSpace& label_space,
      std::span<const Mask> child_rows,
      const std::vector<Mask>& parent_invariant)
      : physical_words_(all_codewords(child_rows)),
        stabilizer_words_(all_codewords(label_space.stabilizer_basis)),
        profiles_(physical_words_.size()),
        parent_invariant_(parent_invariant),
        profile_size_(label_space.points.size() + 1) {
    // Functional coefficients refer to the caller's basis, not its RREF.
    for (auto& word : physical_words_) {
      if ((word >> label_space.quotient_dimension()) != 0) {
        throw std::invalid_argument("quotient row lies outside label space");
      }
      word = linear_combination(word, label_space.quotient_basis);
    }
    const auto header = static_cast<Mask>(label_space.points.size()) |
        (static_cast<Mask>(label_space.stabilizer_basis.size()) << 8U) |
        (static_cast<Mask>(child_rows.size() - 1) << 16U);
    if (parent_invariant_.size() !=
            1 + (physical_words_.size() / 2 - 1) * profile_size_ ||
        parent_invariant_.front() != header) {
      throw std::logic_error("hyperplane invariant has inconsistent dimensions");
    }
  }

  [[nodiscard]] bool smaller(Mask functional) {
    std::vector<Mask> indices;
    indices.reserve(physical_words_.size() / 2 - 1);
    for (Mask coefficients = 1; coefficients < physical_words_.size();
         ++coefficients) {
      if ((std::popcount(coefficients & functional) & 1U) != 0) continue;
      auto& profile = profiles_[coefficients];
      // Different hyperplanes share cosets of the fixed stabiliser space.
      if (profile.empty()) {
        profile.assign(profile_size_, 0);
        for (const auto stabilizer : stabilizer_words_) {
          ++profile[std::popcount(physical_words_[coefficients] ^ stabilizer)];
        }
        for (auto& count : profile) count = stabilizer_words_.size() - count;
      }
      indices.push_back(coefficients);
    }
    std::sort(indices.begin(), indices.end(), [&](Mask left, Mask right) {
      return profiles_[left] < profiles_[right];
    });
    std::size_t offset = 1;
    for (const auto index : indices) {
      for (const auto count : profiles_[index]) {
        const auto previous = parent_invariant_[offset++];
        if (count != previous) return count < previous;
      }
    }
    return false;
  }

 private:
  std::vector<Mask> physical_words_;
  std::vector<Mask> stabilizer_words_;
  std::vector<std::vector<Mask>> profiles_;
  const std::vector<Mask>& parent_invariant_;
  std::size_t profile_size_;
};

[[nodiscard]] bool has_smaller_hyperplane_invariant(
    const LabelSpace& label_space,
    std::span<const Mask> child_rows,
    const std::vector<Mask>& parent_invariant,
    bool nondegenerate_only,
    std::uint32_t nondegenerate_seed_dimension = 0) {
  const auto dimension = static_cast<std::uint32_t>(child_rows.size());
  if (dimension <= 1) return false;
  HyperplaneInvariantCache invariant_cache(label_space, child_rows, parent_invariant);
  const auto functional_count = Mask{1} << dimension;
  for (Mask functional = 1; functional < functional_count; ++functional) {
    const std::array<Mask, 1> constraint{functional};
    const auto kernel = nullspace_basis(constraint, dimension);
    std::vector<Mask> hyperplane;
    hyperplane.reserve(kernel.size());
    for (const auto coefficients : kernel) {
      hyperplane.push_back(linear_combination(coefficients, child_rows));
    }
    if (nondegenerate_only &&
        !tensor_radical_basis(
             physical_label_rows(label_space, hyperplane)).empty()) {
      continue;
    }
    if (!invariant_cache.smaller(functional)) {
      continue;
    }
    // Intermediate parents must retain the seed used to reach this level.
    if (nondegenerate_seed_dimension != 0 &&
        !has_nondegenerate_tensor_restriction(
            dimension - 1,
            cubic_tensor_word(physical_label_rows(label_space, hyperplane)),
            nondegenerate_seed_dimension)) {
      continue;
    }
    return true;
  }
  return false;
}

enum class Q5ParentDomain : std::uint8_t {
  unrestricted,
  q3_chain_cover,
  q4_hitting_set,
};

[[nodiscard]] const std::array<bool, 128>& q5_chain_q3_targets() {
  static const std::array<bool, 128> membership = [] {
    std::array<bool, 128> result{};
    const auto d3_03 =
        canonicalize_cubic_tensor_direct(3, 0x18ULL).canonical_key_words;
    const auto d3_05 =
        canonicalize_cubic_tensor_direct(3, 0x41ULL).canonical_key_words;
    for (std::uint64_t signature = 0; signature < result.size(); ++signature) {
      const auto key =
          canonicalize_cubic_tensor_direct(3, signature).canonical_key_words;
      result[signature] = key == d3_03 || key == d3_05;
    }
    return result;
  }();
  return membership;
}

[[nodiscard]] const std::array<bool, 1U << 14U>&
q5_chain_q4_contains_target() {
  static const std::array<bool, 1U << 14U> membership = [] {
    std::array<bool, 1U << 14U> result{};
    const auto& q3_targets = q5_chain_q3_targets();
    for (std::uint64_t signature = 0; signature < result.size(); ++signature) {
      for (Mask functional = 1; functional < (Mask{1} << 4U); ++functional) {
        const std::array<Mask, 1> constraint{functional};
        const auto kernel = nullspace_basis(constraint, 4);
        std::vector<std::uint32_t> basis;
        basis.reserve(kernel.size());
        for (const auto vector : kernel) {
          basis.push_back(static_cast<std::uint32_t>(vector));
        }
        const auto restricted =
            restrict_cubic_tensor_signature(4, signature, basis);
        if (q3_targets[restricted]) {
          result[signature] = true;
          break;
        }
      }
    }
    return result;
  }();
  return membership;
}

[[nodiscard]] bool q5_chain_parent_is_eligible(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows) {
  const auto dimension = quotient_rows.size();
  if (dimension != 3 && dimension != 4) {
    throw std::logic_error("q5 chain parent has unexpected dimension");
  }
  const auto label_rows = physical_label_rows(label_space, quotient_rows);
  const auto signature = cubic_tensor_word(label_rows);
  if (dimension == 3) return q5_chain_q3_targets()[signature];
  return cubic_tensor_radical_basis_from_signature(4, signature).empty() &&
         q5_chain_q4_contains_target()[signature];
}

[[nodiscard]] bool has_smaller_q5_chain_hyperplane_invariant(
    const LabelSpace& label_space,
    std::span<const Mask> child_rows,
    const std::vector<Mask>& parent_invariant) {
  const auto dimension = static_cast<std::uint32_t>(child_rows.size());
  if (dimension != 4 && dimension != 5) {
    throw std::logic_error("q5 chain child has unexpected dimension");
  }
  const auto functional_count = Mask{1} << dimension;
  for (Mask functional = 1; functional < functional_count; ++functional) {
    const std::array<Mask, 1> constraint{functional};
    const auto kernel = nullspace_basis(constraint, dimension);
    std::vector<Mask> hyperplane;
    hyperplane.reserve(kernel.size());
    for (const auto coefficients : kernel) {
      hyperplane.push_back(linear_combination(coefficients, child_rows));
    }
    if (!q5_chain_parent_is_eligible(label_space, hyperplane)) continue;
    if (marked_code_invariant(label_space, hyperplane) < parent_invariant) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] const std::array<bool, 1U << 14U>&
q5_q4_hitting_set_targets() {
  static const std::array<bool, 1U << 14U> membership = [] {
    constexpr std::array<std::uint64_t, 7> representatives{
        0x000fULL,
        0x0070ULL,
        0x0071ULL,
        0x00bfULL,
        0x00c0ULL,
        0x00f6ULL,
        0x0408ULL,
    };
    std::array<bool, 1U << 14U> result{};
    std::vector<std::uint64_t> queue;
    for (const auto representative : representatives) {
      if (result[representative]) {
        throw std::logic_error("q5 q4 hitting-set targets are not distinct");
      }
      queue.clear();
      queue.push_back(representative);
      result[representative] = true;
      for (std::size_t head = 0; head < queue.size(); ++head) {
        const auto signature = queue[head];
        for (std::uint32_t target = 0; target < 4; ++target) {
          for (std::uint32_t source = 0; source < 4; ++source) {
            if (target == source) continue;
            std::array<std::uint32_t, 4> basis{
                1U, 2U, 4U, 8U};
            basis[target] ^= basis[source];
            const auto image =
                transform_cubic_tensor_signature(4, signature, basis);
            if (!result[image]) {
              result[image] = true;
              queue.push_back(image);
            }
          }
        }
      }
    }
    return result;
  }();
  return membership;
}

[[nodiscard]] bool has_smaller_q5_q4_hitting_hyperplane_invariant(
    const LabelSpace& label_space,
    std::span<const Mask> child_rows,
    const std::vector<Mask>& parent_invariant) {
  const auto dimension = static_cast<std::uint32_t>(child_rows.size());
  if (dimension != 5) {
    throw std::logic_error("q5 q4-hitting child has unexpected dimension");
  }
  const auto& targets = q5_q4_hitting_set_targets();
  const auto functional_count = Mask{1} << dimension;
  for (Mask functional = 1; functional < functional_count; ++functional) {
    const std::array<Mask, 1> constraint{functional};
    const auto kernel = nullspace_basis(constraint, dimension);
    std::vector<Mask> hyperplane;
    hyperplane.reserve(kernel.size());
    for (const auto coefficients : kernel) {
      hyperplane.push_back(linear_combination(coefficients, child_rows));
    }
    const auto signature =
        cubic_tensor_word(physical_label_rows(label_space, hyperplane));
    if (!targets[signature]) continue;
    if (marked_code_invariant(label_space, hyperplane) < parent_invariant) {
      return true;
    }
  }
  return false;
}

#ifdef UTSP_HAVE_BLISS

[[nodiscard]] std::vector<Mask> canonical_key(
    std::span<const Mask> stabilizer_basis,
    std::span<const Mask> full_basis,
    std::span<const unsigned int> canonical_label,
    std::size_t coordinate_count,
    std::size_t stabilizer_dimension,
    std::size_t logical_dimension) {
  std::vector<std::pair<unsigned int, std::uint32_t>> ordered_coordinates;
  ordered_coordinates.reserve(coordinate_count);
  for (std::uint32_t coordinate = 0; coordinate < coordinate_count;
       ++coordinate) {
    ordered_coordinates.emplace_back(canonical_label[coordinate], coordinate);
  }
  std::sort(ordered_coordinates.begin(), ordered_coordinates.end());
  std::vector<std::uint32_t> coordinate_images(coordinate_count);
  for (std::uint32_t image = 0; image < ordered_coordinates.size(); ++image) {
    coordinate_images[ordered_coordinates[image].second] = image;
  }

  std::vector<Mask> canonical_stabilizer;
  canonical_stabilizer.reserve(stabilizer_basis.size());
  for (const auto row : stabilizer_basis) {
    canonical_stabilizer.push_back(permute_row(row, coordinate_images));
  }
  canonical_stabilizer = rref_basis(canonical_stabilizer);
  std::vector<Mask> canonical_full;
  canonical_full.reserve(full_basis.size());
  for (const auto row : full_basis) {
    canonical_full.push_back(permute_row(row, coordinate_images));
  }
  canonical_full = rref_basis(canonical_full);

  std::vector<Mask> key;
  key.reserve(
      1 + canonical_stabilizer.size() + canonical_full.size());
  key.push_back(
      static_cast<Mask>(coordinate_count) |
      (static_cast<Mask>(stabilizer_dimension) << 8U) |
      (static_cast<Mask>(logical_dimension) << 16U));
  key.insert(
      key.end(), canonical_stabilizer.begin(), canonical_stabilizer.end());
  key.insert(
      key.end(),
      canonical_full.begin(),
      canonical_full.end());
  return key;
}

#endif

struct ExtensionQuotient {
  std::vector<Mask> lifts;
  std::vector<std::vector<Mask>> generator_images;
};

[[nodiscard]] ExtensionQuotient build_extension_quotient(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows,
    std::span<const std::vector<Mask>> automorphism_generators) {
  const auto rows = rref_basis(quotient_rows);
  const auto width = label_space.quotient_dimension();
  std::vector<Mask> constraints;
  constraints.reserve(
      rows.size() * label_space.bilinear_form_rows.size());
  for (const auto row : rows) {
    for (const auto& form : label_space.bilinear_form_rows) {
      constraints.push_back(linear_combination(row, form));
    }
  }
  const auto orthogonal_basis = nullspace_basis(constraints, width);

  ExtensionQuotient result;
  std::vector<Mask> span_basis = rows;
  for (const auto vector : orthogonal_basis) {
    const auto residual = reduce_vector(vector, span_basis);
    if (residual == 0) continue;
    result.lifts.push_back(residual);
    span_basis.push_back(residual);
    span_basis = rref_basis(span_basis);
  }
  std::vector<Mask> coordinate_basis = rows;
  coordinate_basis.insert(
      coordinate_basis.end(), result.lifts.begin(), result.lifts.end());
  result.generator_images.reserve(automorphism_generators.size());
  for (const auto& generator : automorphism_generators) {
    if (generator.size() != width) {
      throw std::logic_error("quotient automorphism has the wrong dimension");
    }
    std::vector<Mask> induced;
    induced.reserve(result.lifts.size());
    for (const auto lift : result.lifts) {
      const auto image = linear_combination(lift, generator);
      const auto coordinates = coordinates_in_basis(image, coordinate_basis);
      if (!coordinates) {
        throw std::logic_error(
            "subspace automorphism did not preserve its common orthogonal");
      }
      induced.push_back(*coordinates >> rows.size());
    }
    if (gf2_rank(induced) != result.lifts.size()) {
      throw std::logic_error("extension-quotient automorphism is singular");
    }
    result.generator_images.push_back(std::move(induced));
  }
  return result;
}

[[nodiscard]] std::vector<Mask> vector_orbit_representatives(
    std::uint32_t dimension,
    std::span<const std::vector<Mask>> generators) {
  if (dimension >= 31) {
    throw std::invalid_argument("explicit extension orbits require dimension <31");
  }
  const auto count = std::uint32_t{1} << dimension;
  std::vector<bool> seen(count, false);
  std::vector<Mask> queue;
  std::vector<Mask> representatives;
  seen[0] = true;
  for (Mask seed = 1; seed < count; ++seed) {
    if (seen[seed]) continue;
    representatives.push_back(seed);
    queue.clear();
    queue.push_back(seed);
    seen[seed] = true;
    for (std::size_t next = 0; next < queue.size(); ++next) {
      for (const auto& generator : generators) {
        const auto image = linear_combination(queue[next], generator);
        if (seen[image]) continue;
        seen[image] = true;
        queue.push_back(image);
      }
    }
  }
  return representatives;
}

struct VectorHash {
  [[nodiscard]] std::size_t operator()(
      const std::vector<Mask>& values) const noexcept {
    std::size_t state = 1469598103934665603ULL;
    for (const auto value : values) {
      state ^= static_cast<std::size_t>(value);
      state *= 1099511628211ULL;
    }
    return state;
  }
};

struct CanonicalOrbitNode {
  std::vector<Mask> rows;
  MarkedCodeCanonicalForm canonical;
};

// Local orbit keys are valid only for this fixed support. Public keys still
// come from Bliss, after the final level has been reached.
struct FixedSupportOrbitNode {
  std::vector<Mask> rows;
  MarkedCodeCanonicalForm canonical;
  std::vector<Mask> local_key;
  Mask small_stabilizer_mask = 0;
  std::vector<std::uint32_t> stabilizer_indices;
};

class SmallQuotientGroup {
 public:
  [[nodiscard]] static std::optional<SmallQuotientGroup> build(
      const MarkedCodeCanonicalForm& root,
      std::uint32_t width,
      std::size_t maximum_order = 4096) {
    if (maximum_order == 0) return std::nullopt;
    if (maximum_order > 4096) {
      throw std::invalid_argument("explicit quotient group cap exceeds 4096 elements");
    }
    SmallQuotientGroup result;
    std::vector<Mask> identity;
    for (std::uint32_t bit = 0; bit < width; ++bit) {
      identity.push_back(Mask{1} << bit);
    }
    result.elements_.push_back(identity);
    std::set<std::vector<Mask>> seen{identity};
    for (std::size_t index = 0; index < result.elements_.size(); ++index) {
      for (const auto& generator : root.quotient_automorphism_generators) {
        if (generator.size() != width) {
          throw std::logic_error("support quotient generator has the wrong width");
        }
        std::vector<Mask> product;
        product.reserve(width);
        for (const auto image : result.elements_[index]) {
          product.push_back(linear_combination(image, generator));
        }
        if (!seen.insert(product).second) continue;
        if (result.elements_.size() == maximum_order) return std::nullopt;
        result.elements_.push_back(std::move(product));
      }
    }
    if (root.automorphism_group_order == 0 ||
        root.automorphism_group_order % result.elements_.size() != 0) {
      throw std::logic_error("quotient group image order does not divide support group");
    }
    result.kernel_order_ = root.automorphism_group_order / result.elements_.size();
    return result;
  }

  [[nodiscard]] const std::vector<std::vector<Mask>>& elements() const {
    return elements_;
  }

  [[nodiscard]] bool accepts_line_parent(
      Mask parent, std::span<const Mask> child,
      std::span<const Mask> canonical_child) const {
    if (child.size() != 2 || canonical_child.size() != 2 || parent == 0) {
      throw std::invalid_argument("line-parent selector requires an incident line and two-space");
    }
    if (parent != child[0] && parent != child[1] && parent != (child[0] ^ child[1])) {
      throw std::invalid_argument("line parent is not contained in the two-space");
    }
    const auto least_line = std::min({canonical_child[0], canonical_child[1],
                                     canonical_child[0] ^ canonical_child[1]});
    for (const auto& element : elements_) {
      if (linear_combination(parent, element) != least_line) continue;
      const std::array<Mask, 2> image{linear_combination(child[0], element),
                                    linear_combination(child[1], element)};
      const auto reduced = rref_basis(image);
      if (std::equal(reduced.begin(), reduced.end(), canonical_child.begin(), canonical_child.end())) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] std::vector<std::vector<Mask>> stabilizer_generators(
      Mask membership) const {
    if (elements_.size() > 64 || (membership & 1U) == 0 ||
        (elements_.size() < 64 && (membership >> elements_.size()) != 0)) {
      throw std::logic_error("invalid small quotient stabilizer membership");
    }
    std::vector<std::vector<Mask>> generators;
    generators.reserve(std::popcount(membership) - 1);
    membership &= ~Mask{1};
    while (membership != 0) {
      const auto index = static_cast<std::size_t>(std::countr_zero(membership));
      generators.push_back(elements_[index]);
      membership &= membership - 1;
    }
    return generators;
  }

  [[nodiscard]] std::vector<std::vector<Mask>> stabilizer_generators(
      const FixedSupportOrbitNode& node) const {
    if (elements_.size() <= 64) return stabilizer_generators(node.small_stabilizer_mask);
    if (node.stabilizer_indices.empty() || node.stabilizer_indices.front() != 0) {
      throw std::logic_error("explicit quotient stabilizer is missing the identity");
    }
    std::vector<std::vector<Mask>> generators;
    generators.reserve(node.stabilizer_indices.size() - 1);
    for (std::size_t index = 1; index < node.stabilizer_indices.size(); ++index) {
      generators.push_back(elements_.at(node.stabilizer_indices[index]));
    }
    return generators;
  }

  [[nodiscard]] FixedSupportOrbitNode classify(
      std::span<const Mask> rows) const {
    FixedSupportOrbitNode result{rref_basis(rows), {}, {}, 0, {}};
    result.local_key = result.rows;
    std::uint64_t stabilizer_order = 0;
    for (std::size_t index = 0; index < elements_.size(); ++index) {
      std::vector<Mask> image;
      image.reserve(result.rows.size());
      for (const auto row : result.rows) {
        image.push_back(linear_combination(row, elements_[index]));
      }
      image = rref_basis(image);
      if (image < result.local_key) result.local_key = image;
      if (image != result.rows) continue;
      ++stabilizer_order;
      if (elements_.size() <= 64) {
        result.small_stabilizer_mask |= Mask{1} << index;
      } else {
        result.stabilizer_indices.push_back(static_cast<std::uint32_t>(index));
      }
    }
    result.canonical.automorphism_group_order = kernel_order_ * stabilizer_order;
    return result;
  }

 private:
  std::vector<std::vector<Mask>> elements_;
  std::uint64_t kernel_order_ = 0;
};

}  // namespace

bool marked_code_canonicalizer_available() {
#ifdef UTSP_HAVE_BLISS
  return true;
#else
  return false;
#endif
}

MarkedCodeCanonicalForm canonicalize_marked_code(
    const LabelSpace& label_space,
    std::span<const Mask> quotient_rows) {
#ifndef UTSP_HAVE_BLISS
  (void)label_space;
  (void)quotient_rows;
  throw std::runtime_error(
      "marked-code canonicalization requires a Bliss-enabled build");
#else
  const auto canonical_rows = rref_basis(quotient_rows);
  if (canonical_rows.size() != quotient_rows.size()) {
    throw std::invalid_argument("logical quotient rows are linearly dependent");
  }
  const auto logical_rows = physical_label_rows(label_space, canonical_rows);
  std::vector<Mask> full_basis = label_space.stabilizer_basis;
  full_basis.insert(full_basis.end(), logical_rows.begin(), logical_rows.end());
  if (gf2_rank(full_basis) != full_basis.size()) {
    throw std::logic_error("marked code basis unexpectedly lost rank");
  }
  if (full_basis.size() >= 24) {
    throw std::invalid_argument(
        "marked-code graph is limited to code dimension below 24");
  }

  const auto codewords = all_codewords(full_basis);
  const auto stabilizer_count =
      std::size_t{1} << label_space.stabilizer_basis.size();
  const auto stabilizer_words =
      std::span<const Mask>(codewords).first(stabilizer_count);
  const auto nonstabilizer_words =
      std::span<const Mask>(codewords).subspan(stabilizer_count);

  bliss::Graph graph;
  const auto coordinate_count = label_space.points.size();
  for (std::size_t coordinate = 0; coordinate < coordinate_count; ++coordinate) {
    graph.add_vertex(0);
  }
  std::vector<unsigned int> stabilizer_vertices;
  stabilizer_vertices.reserve(stabilizer_words.size());
  for (std::size_t index = 0; index < stabilizer_words.size(); ++index) {
    stabilizer_vertices.push_back(graph.add_vertex(1));
  }
  std::vector<unsigned int> nonstabilizer_vertices;
  nonstabilizer_vertices.reserve(nonstabilizer_words.size());
  for (std::size_t index = 0; index < nonstabilizer_words.size(); ++index) {
    nonstabilizer_vertices.push_back(graph.add_vertex(2));
  }
  auto add_incidence_edges = [&](std::span<const Mask> words,
                                 std::span<const unsigned int> vertices) {
    for (std::size_t index = 0; index < words.size(); ++index) {
      Mask word = words[index];
      while (word != 0) {
        const auto coordinate = static_cast<unsigned int>(
            std::countr_zero(word));
        graph.add_edge(vertices[index], coordinate);
        word &= word - 1;
      }
    }
  };
  add_incidence_edges(stabilizer_words, stabilizer_vertices);
  add_incidence_edges(nonstabilizer_words, nonstabilizer_vertices);
  graph.set_splitting_heuristic(bliss::Graph::shs_fsm);

  MarkedCodeCanonicalForm result;
  std::set<std::vector<Mask>> generator_set;
  bliss::Stats statistics;
  const auto* labeling = graph.canonical_form(
      statistics,
      [&](unsigned int vertex_count, const unsigned int* automorphism) {
        if (vertex_count < coordinate_count) {
          throw std::logic_error("Bliss automorphism has the wrong degree");
        }
        std::vector<std::uint32_t> coordinate_images(coordinate_count);
        bool identity = true;
        for (std::uint32_t coordinate = 0; coordinate < coordinate_count;
             ++coordinate) {
          if (automorphism[coordinate] >= coordinate_count) {
            throw std::logic_error(
                "Bliss automorphism did not preserve coordinate color");
          }
          coordinate_images[coordinate] = automorphism[coordinate];
          identity &= coordinate_images[coordinate] == coordinate;
        }
        if (identity) return;
        generator_set.insert(induce_quotient_map(
            label_space, coordinate_images));
      });
  if (labeling == nullptr) {
    throw std::runtime_error("Bliss did not return a canonical labeling");
  }
  result.key_words = canonical_key(
      label_space.stabilizer_basis,
      full_basis,
      std::span<const unsigned int>(labeling, graph.get_nof_vertices()),
      coordinate_count,
      label_space.stabilizer_basis.size(),
      canonical_rows.size());
  result.quotient_automorphism_generators.assign(
      generator_set.begin(), generator_set.end());
  const auto approximate_order = statistics.get_group_size_approx();
  if (approximate_order < 1.0L ||
      approximate_order >
          static_cast<long double>(std::numeric_limits<std::uint64_t>::max())) {
    throw std::overflow_error("marked-code automorphism group order overflow");
  }
  result.automorphism_group_order = static_cast<std::uint64_t>(
      std::llround(approximate_order));
  result.canonical_search_nodes = statistics.get_nof_nodes();
  return result;
#endif
}

[[nodiscard]] IsotropicOrbitLevel enumerate_orbits_impl(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t minimum_nondegenerate_dimension,
    bool primitive_at_final_dimension,
    const TensorRestrictionProfile* target_profile,
    std::uint32_t target_filter_from_dimension,
    std::uint32_t exact_nondegenerate_dimension,
    std::uint32_t requested_workers,
    Q5ParentDomain q5_parent_domain = Q5ParentDomain::unrestricted,
    std::set<std::vector<Mask>>* final_tensor_keys = nullptr,
    bool use_small_quotient_group = true,
    bool use_fused_pair_seed = true) {
  if (!marked_code_canonicalizer_available()) {
    throw std::runtime_error(
        "isotropic orbit enumeration requires a Bliss-enabled build");
  }
  if (dimension == 0 || dimension > label_space.quotient_dimension()) {
    throw std::invalid_argument("logical dimension is outside label quotient");
  }

  IsotropicOrbitLevel result;
  result.dimension = dimension;
  result.used_complete_primitive_target_recognizer =
      target_profile != nullptr &&
      target_profile->covers_all_primitive_orbits;
  std::vector<FixedSupportOrbitNode> current;
  current.push_back(FixedSupportOrbitNode{
      {}, canonicalize_marked_code(label_space, {}), {}, 0, {}});
  result.support_automorphism_group_order =
      current.front().canonical.automorphism_group_order;
  result.canonical_search_nodes +=
      current.front().canonical.canonical_search_nodes;
  const auto small_group = use_small_quotient_group
      ? SmallQuotientGroup::build(current.front().canonical,
                                 label_space.quotient_dimension())
      : std::nullopt;
  if (small_group) {
    auto root_stabilizer = small_group->classify({});
    current.front().small_stabilizer_mask = root_stabilizer.small_stabilizer_mask;
    current.front().stabilizer_indices = std::move(root_stabilizer.stabilizer_indices);
  }
  if (std::getenv("UTSP_ORBIT_PROGRESS") != nullptr) {
    std::cerr << "explicit_quotient_image_order="
              << (small_group ? small_group->elements().size() : 0) << '\n';
  }
  for (std::uint32_t level = 0; level < dimension; ++level) {
    if (std::getenv("UTSP_ORBIT_PROGRESS") != nullptr) {
      std::cerr << "orbit_level=" << level << " parents=" << current.size()
                << " quotient_dimension=" << label_space.quotient_dimension() << '\n';
    }
    const bool fuse_pairs = use_fused_pair_seed && small_group && level == 1 &&
        dimension >= 3 && exact_nondegenerate_dimension == 3 &&
        minimum_nondegenerate_dimension == 0 && target_profile != nullptr &&
        target_filter_from_dimension == 1;
    result.parent_orbits += current.size();
    struct WorkerResult {
      std::uint64_t streamed_pair_parents = 0;
      std::uint64_t extension_vector_orbits = 0;
      std::uint64_t filtered_degenerate_children = 0;
      std::uint64_t filtered_nonprimitive_children = 0;
      std::uint64_t filtered_target_children = 0;
      std::uint64_t canonical_parent_rejections = 0;
      std::uint64_t duplicate_children = 0;
      std::uint64_t canonical_search_nodes = 0;
      std::set<std::vector<Mask>> final_tensor_keys;
      std::vector<FixedSupportOrbitNode> children;
      std::exception_ptr error;
    };
    const auto worker_count = std::max<std::uint32_t>(
        1, std::min<std::uint32_t>(requested_workers, current.size()));
    std::atomic<std::size_t> next_parent{0};
    std::vector<WorkerResult> worker_results(worker_count);
    std::vector<std::thread> threads;
    threads.reserve(worker_count);
    for (std::uint32_t worker = 0; worker < worker_count; ++worker) {
      threads.emplace_back([&, worker] {
        auto& local = worker_results[worker];
        std::unordered_map<std::vector<Mask>, std::size_t, VectorHash> seen;
        try {
          std::function<void(const FixedSupportOrbitNode&, std::uint32_t)> process_parent;
          process_parent = [&](const FixedSupportOrbitNode& parent, std::uint32_t parent_level) {
            const auto level = parent_level;
            const auto child_dimension = level + 1;
            const bool extending_past_target = target_profile != nullptr &&
                dimension > target_profile->target_dimension;
            const bool filter_child =
                (minimum_nondegenerate_dimension != 0 && child_dimension >= minimum_nondegenerate_dimension) ||
                child_dimension == exact_nondegenerate_dimension ||
                (extending_past_target && child_dimension > target_profile->target_dimension);
            const bool filter_target = target_profile != nullptr &&
                child_dimension >= target_filter_from_dimension &&
                child_dimension <= target_profile->target_dimension;
            const bool filtered_parent_domain =
                (minimum_nondegenerate_dimension != 0 && level >= minimum_nondegenerate_dimension) ||
                level == exact_nondegenerate_dimension ||
                (extending_past_target && level >= target_profile->target_dimension);
            if (fuse_pairs && level == 2) ++local.streamed_pair_parents;
            const auto parent_invariant = marked_code_invariant(
                label_space, parent.rows);
            const auto parent_label_rows =
                physical_label_rows(label_space, parent.rows);
            const CubicExtensionSignature extension_signature(
                parent_label_rows);
            // The root already carries a short generating set from Bliss.
            const bool explicit_stabilizer = small_group && level != 0;
            const auto small_generators = explicit_stabilizer
                ? small_group->stabilizer_generators(parent)
                : std::vector<std::vector<Mask>>{};
            const auto extension = build_extension_quotient(
                label_space,
                parent.rows,
                explicit_stabilizer ? small_generators
                            : parent.canonical.quotient_automorphism_generators);
            std::vector<Mask> extension_label_lifts;
            extension_label_lifts.reserve(extension.lifts.size());
            for (const auto lift : extension.lifts) {
              extension_label_lifts.push_back(linear_combination(
                  lift, label_space.quotient_basis));
            }
            const auto vector_representatives = vector_orbit_representatives(
                static_cast<std::uint32_t>(extension.lifts.size()),
                extension.generator_images);
            local.extension_vector_orbits += vector_representatives.size();
            for (const auto representative : vector_representatives) {
              std::vector<Mask> label_rows;
              bool label_rows_ready = false;
              Mask extension_label = 0;
              std::optional<std::uint64_t> packed_signature;
              if (filter_child || filter_target) {
                extension_label = linear_combination(
                    representative, extension_label_lifts);
                if (extension_signature.available()) {
                  packed_signature = extension_signature.extend(
                      extension_label);
                } else {
                  label_rows = parent_label_rows;
                  label_rows.push_back(extension_label);
                  label_rows_ready = true;
                }
              }
              const auto ensure_label_rows = [&]()
                  -> const std::vector<Mask>& {
                if (!label_rows_ready) {
                  label_rows = parent_label_rows;
                  label_rows.push_back(extension_label);
                  label_rows_ready = true;
                }
                return label_rows;
              };
              if (filter_child) {
                const auto degenerate = packed_signature
                                            ? !cubic_tensor_radical_basis_from_signature(
                                                   child_dimension,
                                                   *packed_signature)
                                                   .empty()
                                            : !tensor_radical_basis(
                                                   ensure_label_rows())
                                                   .empty();
                if (degenerate) {
                  ++local.filtered_degenerate_children;
                  continue;
                }
                if (primitive_at_final_dimension &&
                    child_dimension == dimension &&
                    (packed_signature
                         ? nondegenerate_tensor_hyperplane_count_from_signature(
                               child_dimension, *packed_signature)
                         : nondegenerate_tensor_hyperplane_count(
                               ensure_label_rows())) !=
                        0) {
                  ++local.filtered_nonprimitive_children;
                  continue;
                }
              }
              if (filter_target) {
                bool allowed = false;
                if (child_dimension == target_profile->target_dimension &&
                    target_profile->covers_all_primitive_orbits) {
                  allowed = is_primitive_exceptional_tensor(
                      ensure_label_rows());
                } else if (
                    target_profile->uses_q7_primitive_family_recognizer) {
                  allowed = target_profile->allows_standard_signature(
                      child_dimension,
                      packed_signature ? *packed_signature
                                       : cubic_tensor_word(label_rows));
                } else {
                  const auto tensor = packed_signature
                                          ? canonicalize_cubic_tensor_direct(
                                                child_dimension,
                                                *packed_signature)
                                          : canonicalize_cubic_tensor_from_labels(
                                                ensure_label_rows());
                  allowed = target_profile->allows(
                      child_dimension, tensor.canonical_key_words);
                }
                if (!allowed) {
                  ++local.filtered_target_children;
                  continue;
                }
              }
              if (final_tensor_keys != nullptr &&
                  child_dimension == dimension) {
                const auto tensor = packed_signature
                                        ? canonicalize_cubic_tensor_direct(
                                              child_dimension,
                                              *packed_signature)
                                        : canonicalize_cubic_tensor_from_labels(
                                              ensure_label_rows());
                local.final_tensor_keys.insert(tensor.canonical_key_words);
                continue;
              }
              const auto extension_row =
                  linear_combination(representative, extension.lifts);
              auto child = parent.rows;
              child.push_back(extension_row);
              bool has_smaller_parent = false;
              const bool use_q5_parent_domain =
                  q5_parent_domain != Q5ParentDomain::unrestricted &&
                  extending_past_target &&
                  level >= target_profile->target_dimension;
              if (fuse_pairs && level == 1) {
                // The exact incidence-orbit selector replaces the old invariant parent rule.
              } else if (!use_q5_parent_domain) {
                has_smaller_parent = has_smaller_hyperplane_invariant(
                    label_space,
                    child,
                    parent_invariant,
                    filtered_parent_domain,
                    target_profile != nullptr &&
                            level > exact_nondegenerate_dimension &&
                            level < target_profile->target_dimension
                        ? exact_nondegenerate_dimension
                        : 0);
              } else if (q5_parent_domain == Q5ParentDomain::q3_chain_cover) {
                has_smaller_parent = has_smaller_q5_chain_hyperplane_invariant(
                    label_space, child, parent_invariant);
              } else {
                has_smaller_parent =
                    has_smaller_q5_q4_hitting_hyperplane_invariant(
                        label_space, child, parent_invariant);
              }
              if (has_smaller_parent) {
                ++local.canonical_parent_rejections;
                continue;
              }
              child = rref_basis(child);
              FixedSupportOrbitNode child_node;
              if (small_group) {
                child_node = small_group->classify(child);
              } else {
                child_node.rows = std::move(child);
                child_node.canonical = canonicalize_marked_code(
                    label_space, child_node.rows);
                child_node.local_key = child_node.canonical.key_words;
              }
              local.canonical_search_nodes +=
                  child_node.canonical.canonical_search_nodes;
              if (fuse_pairs && level == 1) {
                if (!small_group->accepts_line_parent(
                        parent.rows.front(), child_node.rows, child_node.local_key)) {
                  ++local.canonical_parent_rejections;
                  continue;
                }
                process_parent(child_node, 2);
                continue;
              }
              const auto [found, inserted] = seen.emplace(
                  child_node.local_key, local.children.size());
              (void)found;
              if (!inserted) {
                ++local.duplicate_children;
                continue;
              }
              local.children.push_back(std::move(child_node));
            }
          };
          while (true) {
            const auto parent_index = next_parent.fetch_add(1);
            if (parent_index >= current.size()) break;
            process_parent(current[parent_index], level);
          }
        } catch (...) {
          local.error = std::current_exception();
        }
      });
    }
    for (auto& thread : threads) thread.join();

    std::vector<FixedSupportOrbitNode> children;
    for (auto& local : worker_results) {
      if (local.error) std::rethrow_exception(local.error);
      result.parent_orbits += local.streamed_pair_parents;
      result.extension_vector_orbits += local.extension_vector_orbits;
      result.filtered_degenerate_children +=
          local.filtered_degenerate_children;
      result.filtered_nonprimitive_children +=
          local.filtered_nonprimitive_children;
      result.filtered_target_children += local.filtered_target_children;
      result.canonical_parent_rejections +=
          local.canonical_parent_rejections;
      result.duplicate_children += local.duplicate_children;
      result.canonical_search_nodes += local.canonical_search_nodes;
      if (final_tensor_keys != nullptr) {
        final_tensor_keys->insert(
            std::make_move_iterator(local.final_tensor_keys.begin()),
            std::make_move_iterator(local.final_tensor_keys.end()));
      }
      children.insert(
          children.end(),
          std::make_move_iterator(local.children.begin()),
          std::make_move_iterator(local.children.end()));
    }
    std::sort(
        children.begin(), children.end(),
        [](const auto& left, const auto& right) {
          if (left.local_key != right.local_key) {
            return left.local_key < right.local_key;
          }
          return left.rows < right.rows;
        });
    auto retained = children.begin();
    for (auto candidate = children.begin(); candidate != children.end();
         ++candidate) {
      if (retained != children.begin() &&
          (retained - 1)->local_key == candidate->local_key) {
        ++result.duplicate_children;
        continue;
      }
      if (retained != candidate) *retained = std::move(*candidate);
      ++retained;
    }
    children.erase(retained, children.end());
    current = std::move(children);
    if (fuse_pairs) {
      if (std::getenv("UTSP_ORBIT_PROGRESS") != nullptr) {
        std::uint64_t streamed = 0;
        for (const auto& local : worker_results) streamed += local.streamed_pair_parents;
        std::cerr << "streamed_pair_parents=" << streamed << '\n';
      }
      ++level;
    }
  }
  if (small_group) {
    for (auto& node : current) {
      auto canonical = canonicalize_marked_code(label_space, node.rows);
      if (canonical.automorphism_group_order !=
          node.canonical.automorphism_group_order) {
        throw std::logic_error("explicit quotient stabilizer disagrees with Bliss");
      }
      result.canonical_search_nodes += canonical.canonical_search_nodes;
      node.canonical = std::move(canonical);
    }
    std::sort(current.begin(), current.end(), [](const auto& left, const auto& right) {
      return left.canonical.key_words < right.canonical.key_words;
    });
  }
  result.representatives.reserve(current.size());
  result.canonical_keys.reserve(current.size());
  for (auto& node : current) {
    if (node.canonical.automorphism_group_order == 0 ||
        result.support_automorphism_group_order %
                node.canonical.automorphism_group_order !=
            0) {
      throw std::logic_error(
          "marked-code automorphism order does not divide support group");
    }
    const auto orbit_size = result.support_automorphism_group_order /
                            node.canonical.automorphism_group_order;
    if (result.weighted_subspace_count >
        std::numeric_limits<std::uint64_t>::max() - orbit_size) {
      throw std::overflow_error("weighted subspace count overflow");
    }
    result.weighted_subspace_count += orbit_size;
    result.orbit_sizes.push_back(orbit_size);
    result.canonical_keys.push_back(node.canonical.key_words);
    result.representatives.push_back(std::move(node.rows));
  }
  return result;
}

IsotropicOrbitLevel enumerate_isotropic_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  return enumerate_orbits_impl(
      label_space, dimension, 0, false, nullptr, 0, 0, workers);
}

IsotropicOrbitLevel enumerate_hereditary_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  return enumerate_orbits_impl(
      label_space, dimension, 1, false, nullptr, 0, 0, workers);
}

IsotropicOrbitLevel
enumerate_dimension_filtered_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t minimum_nondegenerate_dimension,
    std::uint32_t workers) {
  if (minimum_nondegenerate_dimension == 0 ||
      minimum_nondegenerate_dimension > dimension) {
    throw std::invalid_argument(
        "minimum retained nondegenerate dimension is outside the run");
  }
  return enumerate_orbits_impl(
      label_space,
      dimension,
      minimum_nondegenerate_dimension,
      false,
      nullptr,
      0,
      0,
      workers);
}

IsotropicOrbitLevel enumerate_primitive_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  if (dimension < 2) {
    throw std::invalid_argument("primitive tensor search requires q >= 2");
  }
  return enumerate_orbits_impl(
      label_space, dimension, dimension, true, nullptr, 0, 0, workers);
}

IsotropicOrbitLevel enumerate_target_tensor_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t search_dimension,
    std::uint32_t target_dimension,
    std::span<const std::uint64_t> standard_target_signatures,
    std::uint32_t nondegenerate_seed_dimension,
    std::uint32_t workers) {
  if (search_dimension == 0 || target_dimension == 0 ||
      target_dimension > search_dimension) {
    throw std::invalid_argument(
        "target tensor dimension lies outside the search depth");
  }
  if (nondegenerate_seed_dimension > target_dimension) {
    throw std::invalid_argument(
        "target-tensor nondegenerate seed is deeper than its target");
  }
  if (nondegenerate_seed_dimension != 0) {
    for (const auto signature : standard_target_signatures) {
      if (!has_nondegenerate_tensor_restriction(
              target_dimension,
              signature,
              nondegenerate_seed_dimension)) {
        throw std::invalid_argument(
            "a target tensor has no nondegenerate restriction at the seed dimension");
      }
    }
  }
  // Restriction closure is safe below a required nondegenerate seed as well.
  // Keep the authority-specific q7 recognizer's established start dimension.
  const auto first_filter = target_dimension <= 5 || nondegenerate_seed_dimension == 0
      ? std::uint32_t{1} : nondegenerate_seed_dimension;
  const auto profile = build_tensor_restriction_profile(
      target_dimension,
      standard_target_signatures,
      first_filter,
      workers);
  return enumerate_orbits_impl(
      label_space,
      search_dimension,
      0,
      false,
      &profile,
      first_filter,
      nondegenerate_seed_dimension,
      workers);
}

IsotropicOrbitLevel enumerate_q3_predecessor_target_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 3) {
    throw std::invalid_argument(
        "q3 predecessor target search requires quotient dimension at least 3");
  }

  const std::array<std::uint64_t, 1> d3_03{0x18ULL};
  const std::array<std::uint64_t, 2> primitive_q3{0x40ULL, 0x41ULL};
  std::vector<IsotropicOrbitLevel> branches;
  branches.push_back(enumerate_target_tensor_subspace_orbits(
      label_space, 3, 3, d3_03, 2, workers));
  branches.push_back(enumerate_target_tensor_subspace_orbits(
      label_space, 3, 3, primitive_q3, 3, workers));

  IsotropicOrbitLevel result;
  result.dimension = 3;
  std::map<std::vector<Mask>, std::pair<std::vector<Mask>, std::uint64_t>>
      unique;
  for (auto& branch : branches) {
    if (result.support_automorphism_group_order == 1) {
      result.support_automorphism_group_order =
          branch.support_automorphism_group_order;
    } else if (result.support_automorphism_group_order !=
               branch.support_automorphism_group_order) {
      throw std::logic_error(
          "q3 predecessor branches disagree on support group order");
    }
    result.parent_orbits += branch.parent_orbits;
    result.extension_vector_orbits += branch.extension_vector_orbits;
    result.filtered_degenerate_children +=
        branch.filtered_degenerate_children;
    result.filtered_nonprimitive_children +=
        branch.filtered_nonprimitive_children;
    result.filtered_target_children += branch.filtered_target_children;
    result.canonical_parent_rejections +=
        branch.canonical_parent_rejections;
    result.duplicate_children += branch.duplicate_children;
    result.canonical_search_nodes += branch.canonical_search_nodes;
    result.used_complete_primitive_target_recognizer |=
        branch.used_complete_primitive_target_recognizer;
    if (branch.representatives.size() != branch.canonical_keys.size() ||
        branch.representatives.size() != branch.orbit_sizes.size()) {
      throw std::logic_error("q3 predecessor branch output is misaligned");
    }
    for (std::size_t index = 0; index < branch.representatives.size();
         ++index) {
      const auto [iterator, inserted] = unique.emplace(
          branch.canonical_keys[index],
          std::pair{
              std::move(branch.representatives[index]),
              branch.orbit_sizes[index]});
      if (!inserted) {
        if (iterator->second.second != branch.orbit_sizes[index]) {
          throw std::logic_error(
              "duplicate q3 predecessor orbit has inconsistent size");
        }
        ++result.duplicate_children;
      }
    }
  }
  result.representatives.reserve(unique.size());
  result.canonical_keys.reserve(unique.size());
  result.orbit_sizes.reserve(unique.size());
  for (auto& [key, representative_and_size] : unique) {
    result.weighted_subspace_count += representative_and_size.second;
    result.canonical_keys.push_back(std::move(key));
    result.representatives.push_back(
        std::move(representative_and_size.first));
    result.orbit_sizes.push_back(representative_and_size.second);
  }
  return result;
}

Q3PredecessorExistenceResult find_q3_predecessor_target_subspace(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 3) {
    throw std::invalid_argument(
        "q3 predecessor target search requires quotient dimension at least 3");
  }

  Q3PredecessorExistenceResult result;
  const std::array<std::uint64_t, 1> cs_signature{0x4ULL};
  const auto cs_seeds = enumerate_target_tensor_subspace_orbits(
      label_space, 2, 2, cs_signature, 2, workers);
  result.cs_seed_orbits = cs_seeds.representatives.size();

  static const std::array<bool, 128> is_d3_03_signature = [] {
    std::array<bool, 128> membership{};
    const auto target =
        canonicalize_cubic_tensor_direct(3, 0x18ULL).canonical_key_words;
    for (std::uint64_t signature = 0; signature < membership.size();
         ++signature) {
      membership[signature] =
          canonicalize_cubic_tensor_direct(3, signature)
              .canonical_key_words == target;
    }
    return membership;
  }();

  for (const auto& parent : cs_seeds.representatives) {
    const auto parent_canonical = canonicalize_marked_code(label_space, parent);
    const auto extension = build_extension_quotient(
        label_space,
        parent,
        parent_canonical.quotient_automorphism_generators);
    std::vector<Mask> extension_label_lifts;
    extension_label_lifts.reserve(extension.lifts.size());
    for (const auto lift : extension.lifts) {
      extension_label_lifts.push_back(
          linear_combination(lift, label_space.quotient_basis));
    }
    const auto vector_representatives = vector_orbit_representatives(
        static_cast<std::uint32_t>(extension.lifts.size()),
        extension.generator_images);
    result.extension_vector_orbits += vector_representatives.size();
    const auto parent_label_rows = physical_label_rows(label_space, parent);
    const CubicExtensionSignature signature(parent_label_rows);
    for (const auto representative : vector_representatives) {
      const auto extension_label =
          linear_combination(representative, extension_label_lifts);
      const auto packed = signature.extend(extension_label);
      if (packed >= is_d3_03_signature.size() ||
          !is_d3_03_signature[packed]) {
        continue;
      }
      auto witness = parent;
      witness.push_back(
          linear_combination(representative, extension.lifts));
      result.found = true;
      result.witness_quotient_rows = rref_basis(witness);
      return result;
    }
  }

  const std::array<std::uint64_t, 2> primitive_q3{0x40ULL, 0x41ULL};
  auto primitive = enumerate_target_tensor_subspace_orbits(
      label_space, 3, 3, primitive_q3, 3, workers);
  result.primitive_target_orbits = primitive.representatives.size();
  if (!primitive.representatives.empty()) {
    result.found = true;
    result.witness_quotient_rows =
        std::move(primitive.representatives.front());
  }
  return result;
}

[[nodiscard]] IsotropicOrbitLevel merge_q5_cover_branches(
    std::vector<IsotropicOrbitLevel> branches) {
  IsotropicOrbitLevel result;
  result.dimension = 5;
  std::map<std::vector<Mask>, std::pair<std::vector<Mask>, std::uint64_t>>
      unique;
  for (auto& branch : branches) {
    if (result.support_automorphism_group_order == 1) {
      result.support_automorphism_group_order =
          branch.support_automorphism_group_order;
    } else if (result.support_automorphism_group_order !=
               branch.support_automorphism_group_order) {
      throw std::logic_error(
          "q5 cover branches disagree on support group order");
    }
    result.parent_orbits += branch.parent_orbits;
    result.extension_vector_orbits += branch.extension_vector_orbits;
    result.filtered_degenerate_children +=
        branch.filtered_degenerate_children;
    result.filtered_nonprimitive_children +=
        branch.filtered_nonprimitive_children;
    result.filtered_target_children += branch.filtered_target_children;
    result.canonical_parent_rejections +=
        branch.canonical_parent_rejections;
    result.duplicate_children += branch.duplicate_children;
    result.canonical_search_nodes += branch.canonical_search_nodes;
    result.used_complete_primitive_target_recognizer |=
        branch.used_complete_primitive_target_recognizer;
    if (branch.representatives.size() != branch.canonical_keys.size() ||
        branch.representatives.size() != branch.orbit_sizes.size()) {
      throw std::logic_error("q5 cover branch output is misaligned");
    }
    for (std::size_t index = 0; index < branch.representatives.size();
         ++index) {
      const auto [iterator, inserted] = unique.emplace(
          branch.canonical_keys[index],
          std::pair{
              std::move(branch.representatives[index]),
              branch.orbit_sizes[index]});
      if (!inserted) {
        if (iterator->second.second != branch.orbit_sizes[index]) {
          throw std::logic_error(
              "duplicate q5 cover orbit has inconsistent size");
        }
        ++result.duplicate_children;
      }
    }
  }
  result.representatives.reserve(unique.size());
  result.canonical_keys.reserve(unique.size());
  result.orbit_sizes.reserve(unique.size());
  for (auto& [key, representative_and_size] : unique) {
    result.weighted_subspace_count += representative_and_size.second;
    result.canonical_keys.push_back(std::move(key));
    result.representatives.push_back(
        std::move(representative_and_size.first));
    result.orbit_sizes.push_back(representative_and_size.second);
  }
  return result;
}

IsotropicOrbitLevel enumerate_q5_q3_chain_cover_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 5) {
    throw std::invalid_argument(
        "q5 q3-chain-cover search requires quotient dimension at least 5");
  }

  std::vector<IsotropicOrbitLevel> branches;
  const std::array<std::uint64_t, 1> d3_03{0x18ULL};
  const auto d3_03_profile =
      build_tensor_restriction_profile(3, d3_03, 2, workers);
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &d3_03_profile,
      2,
      2,
      workers,
      Q5ParentDomain::q3_chain_cover));
  const std::array<std::uint64_t, 1> d3_05{0x41ULL};
  const auto d3_05_profile =
      build_tensor_restriction_profile(3, d3_05, 3, workers);
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &d3_05_profile,
      3,
      3,
      workers,
      Q5ParentDomain::q3_chain_cover));
  const std::array<std::uint64_t, 2> primitive_q5{
      0x60000ULL, 0x60001ULL};
  branches.push_back(enumerate_target_tensor_subspace_orbits(
      label_space, 5, 5, primitive_q5, 3, workers));
  return merge_q5_cover_branches(std::move(branches));
}

Q5OutputProfilePresenceResult
enumerate_q5_q3_chain_output_profile_presence(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 5) {
    throw std::invalid_argument(
        "q5 q3-chain profile search requires quotient dimension at least 5");
  }

  std::set<std::vector<Mask>> output_keys;
  std::vector<IsotropicOrbitLevel> branches;
  const std::array<std::uint64_t, 1> d3_03{0x18ULL};
  const auto d3_03_profile =
      build_tensor_restriction_profile(3, d3_03, 2, workers);
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &d3_03_profile,
      2,
      2,
      workers,
      Q5ParentDomain::q3_chain_cover,
      &output_keys));
  const std::array<std::uint64_t, 1> d3_05{0x41ULL};
  const auto d3_05_profile =
      build_tensor_restriction_profile(3, d3_05, 3, workers);
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &d3_05_profile,
      3,
      3,
      workers,
      Q5ParentDomain::q3_chain_cover,
      &output_keys));
  const std::array<std::uint64_t, 2> primitive_q5{
      0x60000ULL, 0x60001ULL};
  const auto primitive_profile =
      build_tensor_restriction_profile(5, primitive_q5, 3, workers);
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &primitive_profile,
      3,
      3,
      workers,
      Q5ParentDomain::unrestricted,
      &output_keys));

  Q5OutputProfilePresenceResult result;
  result.traversal = merge_q5_cover_branches(std::move(branches));
  result.canonical_output_keys.assign(
      std::make_move_iterator(output_keys.begin()),
      std::make_move_iterator(output_keys.end()));
  return result;
}

IsotropicOrbitLevel enumerate_q5_q4_hitting_set_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 5) {
    throw std::invalid_argument(
        "q5 q4-hitting-set search requires quotient dimension at least 5");
  }

  constexpr std::array<std::uint64_t, 7> q4_targets{
      0x000fULL,
      0x0070ULL,
      0x0071ULL,
      0x00bfULL,
      0x00c0ULL,
      0x00f6ULL,
      0x0408ULL,
  };
  const auto q4_profile =
      build_tensor_restriction_profile(4, q4_targets, 1, workers);
  std::vector<IsotropicOrbitLevel> branches;
  branches.push_back(enumerate_orbits_impl(
      label_space,
      5,
      0,
      false,
      &q4_profile,
      1,
      4,
      workers,
      Q5ParentDomain::q4_hitting_set));
  const std::array<std::uint64_t, 2> primitive_q5{
      0x60000ULL, 0x60001ULL};
  branches.push_back(enumerate_target_tensor_subspace_orbits(
      label_space, 5, 5, primitive_q5, 3, workers));
  return merge_q5_cover_branches(std::move(branches));
}

IsotropicOrbitLevel enumerate_complete_nondegenerate_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  if (dimension == 0 || dimension > 8 ||
      dimension > label_space.quotient_dimension()) {
    throw std::invalid_argument(
        "complete nondegenerate orbit search requires q=1,...,8");
  }

  std::vector<IsotropicOrbitLevel> branches;
  branches.push_back(enumerate_hereditary_nondegenerate_subspace_orbits(
      label_space, dimension, workers));
  const auto add_seed = [&](
                            std::uint32_t seed_dimension,
                            std::initializer_list<std::uint64_t> signatures,
                            std::uint32_t nondegenerate_seed) {
    if (seed_dimension > dimension) return;
    const std::vector<std::uint64_t> targets(signatures);
    branches.push_back(enumerate_target_tensor_subspace_orbits(
        label_space,
        dimension,
        seed_dimension,
        targets,
        nondegenerate_seed,
        workers));
  };
  add_seed(2, {0x4ULL}, 2);
  add_seed(3, {0x40ULL, 0x41ULL}, 3);
  add_seed(5, {0x60000ULL, 0x60001ULL}, 3);
  add_seed(7, {0x42010000000ULL, 0x42010000001ULL}, 3);

  IsotropicOrbitLevel result;
  result.dimension = dimension;
  std::map<std::vector<Mask>, std::pair<std::vector<Mask>, std::uint64_t>>
      unique;
  for (auto& branch : branches) {
    if (result.support_automorphism_group_order == 1) {
      result.support_automorphism_group_order =
          branch.support_automorphism_group_order;
    } else if (result.support_automorphism_group_order !=
               branch.support_automorphism_group_order) {
      throw std::logic_error(
          "primitive-seed branches disagree on support group order");
    }
    result.parent_orbits += branch.parent_orbits;
    result.extension_vector_orbits += branch.extension_vector_orbits;
    result.filtered_degenerate_children +=
        branch.filtered_degenerate_children;
    result.filtered_nonprimitive_children +=
        branch.filtered_nonprimitive_children;
    result.filtered_target_children += branch.filtered_target_children;
    result.canonical_parent_rejections +=
        branch.canonical_parent_rejections;
    result.duplicate_children += branch.duplicate_children;
    result.canonical_search_nodes += branch.canonical_search_nodes;
    result.used_complete_primitive_target_recognizer |=
        branch.used_complete_primitive_target_recognizer;
    if (branch.representatives.size() != branch.canonical_keys.size() ||
        branch.representatives.size() != branch.orbit_sizes.size()) {
      throw std::logic_error("primitive-seed branch output is misaligned");
    }
    for (std::size_t index = 0; index < branch.representatives.size();
         ++index) {
      const auto [iterator, inserted] = unique.emplace(
          branch.canonical_keys[index],
          std::pair{
              std::move(branch.representatives[index]),
              branch.orbit_sizes[index]});
      if (!inserted) {
        if (iterator->second.second != branch.orbit_sizes[index]) {
          throw std::logic_error(
              "duplicate primitive-seed orbit has inconsistent size");
        }
        ++result.duplicate_children;
      }
    }
  }
  result.representatives.reserve(unique.size());
  result.canonical_keys.reserve(unique.size());
  result.orbit_sizes.reserve(unique.size());
  for (auto& [key, representative_and_size] : unique) {
    result.weighted_subspace_count += representative_and_size.second;
    result.canonical_keys.push_back(std::move(key));
    result.representatives.push_back(
        std::move(representative_and_size.first));
    result.orbit_sizes.push_back(representative_and_size.second);
  }
  return result;
}

IsotropicOrbitLevel enumerate_nondegenerate_successor_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  if (!marked_code_canonicalizer_available()) {
    throw std::runtime_error(
        "nondegenerate successor search requires a Bliss-enabled build");
  }
  if (dimension < 2 || dimension > 8 ||
      dimension > label_space.quotient_dimension()) {
    throw std::invalid_argument(
        "nondegenerate successor search requires q=2,...,8");
  }

  auto parents = enumerate_complete_nondegenerate_subspace_orbits(
      label_space, dimension - 1, workers);
  IsotropicOrbitLevel result;
  result.dimension = dimension;
  result.parent_orbits =
      parents.parent_orbits + parents.representatives.size();
  result.extension_vector_orbits = parents.extension_vector_orbits;
  result.filtered_degenerate_children = parents.filtered_degenerate_children;
  result.filtered_nonprimitive_children =
      parents.filtered_nonprimitive_children;
  result.filtered_target_children = parents.filtered_target_children;
  result.canonical_parent_rejections = parents.canonical_parent_rejections;
  result.duplicate_children = parents.duplicate_children;
  result.canonical_search_nodes = parents.canonical_search_nodes;
  result.support_automorphism_group_order =
      parents.support_automorphism_group_order;
  result.used_complete_primitive_target_recognizer =
      parents.used_complete_primitive_target_recognizer;

  if (parents.representatives.size() != parents.canonical_keys.size() ||
      parents.representatives.size() != parents.orbit_sizes.size()) {
    throw std::logic_error("nondegenerate successor parents are misaligned");
  }

  std::vector<CanonicalOrbitNode> current;
  current.reserve(parents.representatives.size());
  for (std::size_t index = 0; index < parents.representatives.size(); ++index) {
    auto canonical =
        canonicalize_marked_code(label_space, parents.representatives[index]);
    result.canonical_search_nodes += canonical.canonical_search_nodes;
    if (canonical.key_words != parents.canonical_keys[index]) {
      throw std::logic_error(
          "recanonicalized nondegenerate successor parent changed key");
    }
    if (canonical.automorphism_group_order == 0 ||
        result.support_automorphism_group_order %
                canonical.automorphism_group_order !=
            0 ||
        result.support_automorphism_group_order /
                canonical.automorphism_group_order !=
            parents.orbit_sizes[index]) {
      throw std::logic_error(
          "recanonicalized nondegenerate successor parent changed orbit size");
    }
    current.push_back(CanonicalOrbitNode{
        std::move(parents.representatives[index]), std::move(canonical)});
  }

  struct WorkerResult {
    std::uint64_t extension_vector_orbits = 0;
    std::uint64_t filtered_degenerate_children = 0;
    std::uint64_t canonical_parent_rejections = 0;
    std::uint64_t duplicate_children = 0;
    std::uint64_t canonical_search_nodes = 0;
    std::vector<CanonicalOrbitNode> children;
    std::exception_ptr error;
  };

  const auto worker_count = std::max<std::uint32_t>(
      1, std::min<std::uint32_t>(workers, current.size()));
  std::atomic<std::size_t> next_parent{0};
  std::vector<WorkerResult> worker_results(worker_count);
  std::vector<std::thread> threads;
  threads.reserve(worker_count);
  for (std::uint32_t worker = 0; worker < worker_count; ++worker) {
    threads.emplace_back([&, worker] {
      auto& local = worker_results[worker];
      std::unordered_map<std::vector<Mask>, std::size_t, VectorHash> seen;
      try {
        while (true) {
          const auto parent_index = next_parent.fetch_add(1);
          if (parent_index >= current.size()) break;
          const auto& parent = current[parent_index];
          const auto parent_invariant =
              marked_code_invariant(label_space, parent.rows);
          const auto parent_label_rows =
              physical_label_rows(label_space, parent.rows);
          const CubicExtensionSignature extension_signature(
              parent_label_rows);
          const auto extension = build_extension_quotient(
              label_space,
              parent.rows,
              parent.canonical.quotient_automorphism_generators);
          std::vector<Mask> extension_label_lifts;
          extension_label_lifts.reserve(extension.lifts.size());
          for (const auto lift : extension.lifts) {
            extension_label_lifts.push_back(linear_combination(
                lift, label_space.quotient_basis));
          }
          const auto vector_representatives = vector_orbit_representatives(
              static_cast<std::uint32_t>(extension.lifts.size()),
              extension.generator_images);
          local.extension_vector_orbits += vector_representatives.size();
          for (const auto representative : vector_representatives) {
            const auto extension_label = linear_combination(
                representative, extension_label_lifts);
            bool degenerate = false;
            if (extension_signature.available()) {
              degenerate = !cubic_tensor_radical_basis_from_signature(
                                static_cast<std::uint32_t>(
                                    parent_label_rows.size() + 1),
                                extension_signature.extend(extension_label))
                                .empty();
            } else {
              auto label_rows = parent_label_rows;
              label_rows.push_back(extension_label);
              degenerate = !tensor_radical_basis(label_rows).empty();
            }
            if (degenerate) {
              ++local.filtered_degenerate_children;
              continue;
            }
            const auto extension_row =
                linear_combination(representative, extension.lifts);
            auto child = parent.rows;
            child.push_back(extension_row);
            if (has_smaller_hyperplane_invariant(
                    label_space, child, parent_invariant, true)) {
              ++local.canonical_parent_rejections;
              continue;
            }
            child = rref_basis(child);
            auto child_canonical =
                canonicalize_marked_code(label_space, child);
            local.canonical_search_nodes +=
                child_canonical.canonical_search_nodes;
            const auto [found, inserted] = seen.emplace(
                child_canonical.key_words, local.children.size());
            (void)found;
            if (!inserted) {
              ++local.duplicate_children;
              continue;
            }
            local.children.push_back(CanonicalOrbitNode{
                std::move(child), std::move(child_canonical)});
          }
        }
      } catch (...) {
        local.error = std::current_exception();
      }
    });
  }
  for (auto& thread : threads) thread.join();

  std::vector<CanonicalOrbitNode> children;
  for (auto& local : worker_results) {
    if (local.error) std::rethrow_exception(local.error);
    result.extension_vector_orbits += local.extension_vector_orbits;
    result.filtered_degenerate_children +=
        local.filtered_degenerate_children;
    result.canonical_parent_rejections +=
        local.canonical_parent_rejections;
    result.duplicate_children += local.duplicate_children;
    result.canonical_search_nodes += local.canonical_search_nodes;
    children.insert(
        children.end(),
        std::make_move_iterator(local.children.begin()),
        std::make_move_iterator(local.children.end()));
  }
  std::sort(
      children.begin(), children.end(),
      [](const auto& left, const auto& right) {
        if (left.canonical.key_words != right.canonical.key_words) {
          return left.canonical.key_words < right.canonical.key_words;
        }
        return left.rows < right.rows;
      });
  auto retained = children.begin();
  for (auto candidate = children.begin(); candidate != children.end();
       ++candidate) {
    if (retained != children.begin() &&
        (retained - 1)->canonical.key_words ==
            candidate->canonical.key_words) {
      ++result.duplicate_children;
      continue;
    }
    if (retained != candidate) *retained = std::move(*candidate);
    ++retained;
  }
  children.erase(retained, children.end());

  result.representatives.reserve(children.size());
  result.canonical_keys.reserve(children.size());
  result.orbit_sizes.reserve(children.size());
  for (auto& node : children) {
    if (node.canonical.automorphism_group_order == 0 ||
        result.support_automorphism_group_order %
                node.canonical.automorphism_group_order !=
            0) {
      throw std::logic_error(
          "successor automorphism order does not divide support group");
    }
    const auto orbit_size = result.support_automorphism_group_order /
                            node.canonical.automorphism_group_order;
    if (result.weighted_subspace_count >
        std::numeric_limits<std::uint64_t>::max() - orbit_size) {
      throw std::overflow_error("successor weighted subspace count overflow");
    }
    result.weighted_subspace_count += orbit_size;
    result.orbit_sizes.push_back(orbit_size);
    result.canonical_keys.push_back(std::move(node.canonical.key_words));
    result.representatives.push_back(std::move(node.rows));
  }
  return result;
}

IsotropicOrbitLevel
enumerate_q5_zero_hyperplane_primitive_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (!marked_code_canonicalizer_available()) {
    throw std::runtime_error(
        "zero-hyperplane primitive search requires a Bliss-enabled build");
  }
  if (workers == 0) {
    throw std::invalid_argument(
        "zero-hyperplane primitive worker count must be positive");
  }

  IsotropicOrbitLevel result;
  result.dimension = 5;
  result.used_complete_primitive_target_recognizer = true;
  if (label_space.quotient_dimension() < 5) return result;

  std::vector<CanonicalOrbitNode> current;
  current.push_back(CanonicalOrbitNode{
      {}, canonicalize_marked_code(label_space, {})});
  result.support_automorphism_group_order =
      current.front().canonical.automorphism_group_order;
  result.canonical_search_nodes +=
      current.front().canonical.canonical_search_nodes;

  for (std::uint32_t level = 0; level < 5; ++level) {
    const auto child_dimension = level + 1;
    const auto final_level = child_dimension == 5;
    result.parent_orbits += current.size();
    struct WorkerResult {
      std::uint64_t extension_vector_orbits = 0;
      std::uint64_t filtered_degenerate_children = 0;
      std::uint64_t filtered_nonprimitive_children = 0;
      std::uint64_t filtered_target_children = 0;
      std::uint64_t canonical_parent_rejections = 0;
      std::uint64_t duplicate_children = 0;
      std::uint64_t canonical_search_nodes = 0;
      std::vector<CanonicalOrbitNode> children;
      std::exception_ptr error;
    };
    const auto worker_count = std::max<std::uint32_t>(
        1, std::min<std::uint32_t>(workers, current.size()));
    std::atomic<std::size_t> next_parent{0};
    std::vector<WorkerResult> worker_results(worker_count);
    std::vector<std::thread> threads;
    threads.reserve(worker_count);
    for (std::uint32_t worker = 0; worker < worker_count; ++worker) {
      threads.emplace_back([&, worker] {
        auto& local = worker_results[worker];
        std::unordered_map<std::vector<Mask>, std::size_t, VectorHash> seen;
        try {
          while (true) {
            const auto parent_index = next_parent.fetch_add(1);
            if (parent_index >= current.size()) break;
            const auto& parent = current[parent_index];
            const auto parent_invariant = marked_code_invariant(
                label_space, parent.rows);
            const auto parent_label_rows =
                physical_label_rows(label_space, parent.rows);
            const CubicExtensionSignature extension_signature(
                parent_label_rows);
            const auto extension = build_extension_quotient(
                label_space,
                parent.rows,
                parent.canonical.quotient_automorphism_generators);
            std::vector<Mask> extension_label_lifts;
            extension_label_lifts.reserve(extension.lifts.size());
            for (const auto lift : extension.lifts) {
              extension_label_lifts.push_back(linear_combination(
                  lift, label_space.quotient_basis));
            }
            const auto vector_representatives = vector_orbit_representatives(
                static_cast<std::uint32_t>(extension.lifts.size()),
                extension.generator_images);
            local.extension_vector_orbits += vector_representatives.size();
            for (const auto representative : vector_representatives) {
              const auto extension_label = linear_combination(
                  representative, extension_label_lifts);
              const auto packed_signature =
                  extension_signature.extend(extension_label);
              std::vector<Mask> child;
              std::optional<Mask> extension_row;
              if (!final_level) {
                if (packed_signature != 0) {
                  ++local.filtered_target_children;
                  continue;
                }
                extension_row =
                    linear_combination(representative, extension.lifts);
                child = parent.rows;
                child.push_back(*extension_row);
                if (has_smaller_hyperplane_invariant(
                        label_space, child, parent_invariant, false)) {
                  ++local.canonical_parent_rejections;
                  continue;
                }
              } else {
                if (!cubic_tensor_radical_basis_from_signature(
                         5, packed_signature)
                         .empty()) {
                  ++local.filtered_degenerate_children;
                  continue;
                }
                if (nondegenerate_tensor_hyperplane_count_from_signature(
                        5, packed_signature) != 0) {
                  ++local.filtered_nonprimitive_children;
                  continue;
                }
              }
              if (!extension_row) {
                extension_row =
                    linear_combination(representative, extension.lifts);
              }
              if (child.empty()) {
                child = parent.rows;
                child.push_back(*extension_row);
              }
              child = rref_basis(child);
              auto child_canonical = canonicalize_marked_code(
                  label_space, child);
              local.canonical_search_nodes +=
                  child_canonical.canonical_search_nodes;
              const auto [position, inserted] = seen.emplace(
                  child_canonical.key_words, local.children.size());
              (void)position;
              if (!inserted) {
                ++local.duplicate_children;
                continue;
              }
              local.children.push_back(CanonicalOrbitNode{
                  std::move(child), std::move(child_canonical)});
            }
          }
        } catch (...) {
          local.error = std::current_exception();
        }
      });
    }
    for (auto& thread : threads) thread.join();

    std::vector<CanonicalOrbitNode> children;
    for (auto& local : worker_results) {
      if (local.error) std::rethrow_exception(local.error);
      result.extension_vector_orbits += local.extension_vector_orbits;
      result.filtered_degenerate_children +=
          local.filtered_degenerate_children;
      result.filtered_nonprimitive_children +=
          local.filtered_nonprimitive_children;
      result.filtered_target_children += local.filtered_target_children;
      result.canonical_parent_rejections +=
          local.canonical_parent_rejections;
      result.duplicate_children += local.duplicate_children;
      result.canonical_search_nodes += local.canonical_search_nodes;
      children.insert(
          children.end(),
          std::make_move_iterator(local.children.begin()),
          std::make_move_iterator(local.children.end()));
    }
    std::sort(
        children.begin(), children.end(),
        [](const auto& left, const auto& right) {
          if (left.canonical.key_words != right.canonical.key_words) {
            return left.canonical.key_words < right.canonical.key_words;
          }
          return left.rows < right.rows;
        });
    auto retained = children.begin();
    for (auto candidate = children.begin(); candidate != children.end();
         ++candidate) {
      if (retained != children.begin() &&
          (retained - 1)->canonical.key_words ==
              candidate->canonical.key_words) {
        ++result.duplicate_children;
        continue;
      }
      if (retained != candidate) *retained = std::move(*candidate);
      ++retained;
    }
    children.erase(retained, children.end());
    current = std::move(children);
  }

  result.representatives.reserve(current.size());
  result.canonical_keys.reserve(current.size());
  result.orbit_sizes.reserve(current.size());
  for (auto& node : current) {
    if (node.canonical.automorphism_group_order == 0 ||
        result.support_automorphism_group_order %
                node.canonical.automorphism_group_order !=
            0) {
      throw std::logic_error(
          "zero-hyperplane primitive automorphism order does not divide support group");
    }
    const auto orbit_size = result.support_automorphism_group_order /
                            node.canonical.automorphism_group_order;
    if (result.weighted_subspace_count >
        std::numeric_limits<std::uint64_t>::max() - orbit_size) {
      throw std::overflow_error(
          "zero-hyperplane primitive weighted count overflow");
    }
    result.weighted_subspace_count += orbit_size;
    result.orbit_sizes.push_back(orbit_size);
    result.canonical_keys.push_back(std::move(node.canonical.key_words));
    result.representatives.push_back(std::move(node.rows));
  }
  return result;
}

IsotropicOrbitLevel
enumerate_complete_q5_zero_hyperplane_subspace_orbits(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (label_space.quotient_dimension() < 5) {
    IsotropicOrbitLevel empty;
    empty.dimension = 5;
    return empty;
  }
  const std::array<std::uint64_t, 1> primitive_q2{0x4ULL};
  const std::array<std::uint64_t, 2> primitive_q3{0x40ULL, 0x41ULL};
  std::array<IsotropicOrbitLevel, 4> branches{
      enumerate_hereditary_nondegenerate_subspace_orbits(
          label_space, 5, workers),
      enumerate_target_tensor_subspace_orbits(
          label_space, 5, 2, primitive_q2, 2, workers),
      enumerate_target_tensor_subspace_orbits(
          label_space, 5, 3, primitive_q3, 3, workers),
      enumerate_q5_zero_hyperplane_primitive_subspace_orbits(
          label_space, workers),
  };
  IsotropicOrbitLevel result;
  result.dimension = 5;
  std::map<std::vector<Mask>, std::pair<std::vector<Mask>, std::uint64_t>>
      unique;
  for (auto& branch : branches) {
    if (result.support_automorphism_group_order == 1) {
      result.support_automorphism_group_order =
          branch.support_automorphism_group_order;
    } else if (result.support_automorphism_group_order !=
               branch.support_automorphism_group_order) {
      throw std::logic_error(
          "q5 zero-hyperplane branches disagree on support group order");
    }
    result.parent_orbits += branch.parent_orbits;
    result.extension_vector_orbits += branch.extension_vector_orbits;
    result.filtered_degenerate_children +=
        branch.filtered_degenerate_children;
    result.filtered_nonprimitive_children +=
        branch.filtered_nonprimitive_children;
    result.filtered_target_children += branch.filtered_target_children;
    result.canonical_parent_rejections +=
        branch.canonical_parent_rejections;
    result.duplicate_children += branch.duplicate_children;
    result.canonical_search_nodes += branch.canonical_search_nodes;
    result.used_complete_primitive_target_recognizer |=
        branch.used_complete_primitive_target_recognizer;
    if (branch.representatives.size() != branch.canonical_keys.size() ||
        branch.representatives.size() != branch.orbit_sizes.size()) {
      throw std::logic_error(
          "q5 zero-hyperplane branch output is misaligned");
    }
    for (std::size_t index = 0; index < branch.representatives.size();
         ++index) {
      const auto [position, inserted] = unique.emplace(
          branch.canonical_keys[index],
          std::pair{
              std::move(branch.representatives[index]),
              branch.orbit_sizes[index]});
      if (!inserted) {
        if (position->second.second != branch.orbit_sizes[index]) {
          throw std::logic_error(
              "q5 zero-hyperplane duplicate orbit changed size");
        }
        ++result.duplicate_children;
      }
    }
  }
  result.representatives.reserve(unique.size());
  result.canonical_keys.reserve(unique.size());
  result.orbit_sizes.reserve(unique.size());
  for (auto& [key, representative_and_size] : unique) {
    if (result.weighted_subspace_count >
        std::numeric_limits<std::uint64_t>::max() -
            representative_and_size.second) {
      throw std::overflow_error(
          "q5 zero-hyperplane complete weighted count overflow");
    }
    result.weighted_subspace_count += representative_and_size.second;
    result.canonical_keys.push_back(std::move(key));
    result.representatives.push_back(
        std::move(representative_and_size.first));
    result.orbit_sizes.push_back(representative_and_size.second);
  }
  return result;
}

RadicalStratifiedQ8Result
search_q8_isotropic_subspaces_by_tensor_radical(
    const LabelSpace& label_space,
    std::uint32_t workers) {
  if (!marked_code_canonicalizer_available()) {
    throw std::runtime_error(
        "tensor-radical q=8 search requires a Bliss-enabled build");
  }
  if (workers == 0) {
    throw std::invalid_argument(
        "tensor-radical q=8 worker count must be positive");
  }

  RadicalStratifiedQ8Result result;
  result.seed_orbits_by_intrinsic_dimension.assign(9, 0);
  result.sector_witness_found.assign(9, false);
  if (label_space.quotient_dimension() < 8) return result;

  const auto seed_nodes = [&](std::uint32_t intrinsic_dimension) {
    std::vector<CanonicalOrbitNode> nodes;
    if (intrinsic_dimension == 0) {
      auto canonical = canonicalize_marked_code(label_space, {});
      result.canonical_search_nodes += canonical.canonical_search_nodes;
      nodes.push_back(CanonicalOrbitNode{{}, std::move(canonical)});
    } else {
      auto level = enumerate_complete_nondegenerate_subspace_orbits(
          label_space, intrinsic_dimension, workers);
      result.parent_orbits += level.parent_orbits;
      result.extension_vector_orbits += level.extension_vector_orbits;
      result.duplicate_children += level.duplicate_children;
      result.canonical_search_nodes += level.canonical_search_nodes;
      if (level.representatives.size() != level.canonical_keys.size()) {
        throw std::logic_error(
            "nondegenerate radical-stratification seeds are misaligned");
      }
      nodes.reserve(level.representatives.size());
      for (std::size_t index = 0; index < level.representatives.size();
           ++index) {
        auto canonical = canonicalize_marked_code(
            label_space, level.representatives[index]);
        result.canonical_search_nodes += canonical.canonical_search_nodes;
        if (canonical.key_words != level.canonical_keys[index]) {
          throw std::logic_error(
              "recanonicalized radical-stratification seed changed key");
        }
        nodes.push_back(CanonicalOrbitNode{
            std::move(level.representatives[index]),
            std::move(canonical)});
      }
    }
    result.seed_orbits += nodes.size();
    result.seed_orbits_by_intrinsic_dimension[intrinsic_dimension] =
        nodes.size();
    return nodes;
  };

  const auto sector_has_witness = [&]
      (std::vector<CanonicalOrbitNode> current,
       std::uint32_t intrinsic_dimension) {
    while (!current.empty()) {
      const auto parent_dimension = static_cast<std::uint32_t>(
          current.front().rows.size());
      if (parent_dimension >= 8) return true;
      for (const auto& node : current) {
        if (node.rows.size() != parent_dimension) {
          throw std::logic_error(
              "radical-stratification level mixes dimensions");
        }
      }
      const auto child_dimension = parent_dimension + 1;
      if (child_dimension < intrinsic_dimension) {
        throw std::logic_error(
            "radical-stratification seed is below intrinsic dimension");
      }
      const auto expected_radical_dimension =
          child_dimension - intrinsic_dimension;
      result.parent_orbits += current.size();

      struct WorkerResult {
        std::uint64_t extension_vector_orbits = 0;
        std::uint64_t rejected_radical_dimension = 0;
        std::uint64_t duplicate_children = 0;
        std::uint64_t canonical_search_nodes = 0;
        std::vector<CanonicalOrbitNode> children;
        std::vector<Mask> witness_rows;
        std::exception_ptr error;
      };
      const auto worker_count = std::max<std::uint32_t>(
          1, std::min<std::uint32_t>(workers, current.size()));
      std::atomic<std::size_t> next_parent{0};
      std::atomic<bool> found{false};
      std::vector<WorkerResult> worker_results(worker_count);
      std::vector<std::thread> threads;
      threads.reserve(worker_count);
      for (std::uint32_t worker = 0; worker < worker_count; ++worker) {
        threads.emplace_back([&, worker] {
          auto& local = worker_results[worker];
          std::unordered_map<std::vector<Mask>, std::size_t, VectorHash> seen;
          try {
            while (!found.load(std::memory_order_relaxed)) {
              const auto parent_index = next_parent.fetch_add(1);
              if (parent_index >= current.size()) break;
              const auto& parent = current[parent_index];
              const auto extension = build_extension_quotient(
                  label_space,
                  parent.rows,
                  parent.canonical.quotient_automorphism_generators);
              const auto vector_representatives =
                  vector_orbit_representatives(
                      static_cast<std::uint32_t>(extension.lifts.size()),
                      extension.generator_images);
              local.extension_vector_orbits +=
                  vector_representatives.size();
              for (const auto representative : vector_representatives) {
                if (found.load(std::memory_order_relaxed)) break;
                auto child = parent.rows;
                child.push_back(
                    linear_combination(representative, extension.lifts));
                child = rref_basis(child);
                const auto label_rows =
                    physical_label_rows(label_space, child);
                if (tensor_radical_basis(label_rows).size() !=
                    expected_radical_dimension) {
                  ++local.rejected_radical_dimension;
                  continue;
                }
                if (child_dimension == 8) {
                  bool expected = false;
                  if (found.compare_exchange_strong(
                          expected, true, std::memory_order_relaxed)) {
                    local.witness_rows = std::move(child);
                  }
                  break;
                }
                auto canonical = canonicalize_marked_code(
                    label_space, child);
                local.canonical_search_nodes +=
                    canonical.canonical_search_nodes;
                const auto [position, inserted] = seen.emplace(
                    canonical.key_words, local.children.size());
                (void)position;
                if (!inserted) {
                  ++local.duplicate_children;
                  continue;
                }
                local.children.push_back(CanonicalOrbitNode{
                    std::move(child), std::move(canonical)});
              }
            }
          } catch (...) {
            local.error = std::current_exception();
            found.store(true, std::memory_order_relaxed);
          }
        });
      }
      for (auto& thread : threads) thread.join();

      std::vector<CanonicalOrbitNode> children;
      std::vector<Mask> witness_rows;
      for (auto& local : worker_results) {
        if (local.error) std::rethrow_exception(local.error);
        result.extension_vector_orbits += local.extension_vector_orbits;
        result.rejected_radical_dimension +=
            local.rejected_radical_dimension;
        result.duplicate_children += local.duplicate_children;
        result.canonical_search_nodes += local.canonical_search_nodes;
        if (witness_rows.empty() && !local.witness_rows.empty()) {
          witness_rows = std::move(local.witness_rows);
        }
        children.insert(
            children.end(),
            std::make_move_iterator(local.children.begin()),
            std::make_move_iterator(local.children.end()));
      }
      if (found.load(std::memory_order_relaxed)) {
        if (witness_rows.size() != 8) {
          throw std::logic_error(
              "radical-stratification witness was not retained");
        }
        result.witness_intrinsic_dimension = intrinsic_dimension;
        result.witness_quotient_rows = std::move(witness_rows);
        return true;
      }

      std::sort(
          children.begin(), children.end(),
          [](const auto& left, const auto& right) {
            if (left.canonical.key_words != right.canonical.key_words) {
              return left.canonical.key_words < right.canonical.key_words;
            }
            return left.rows < right.rows;
          });
      auto retained = children.begin();
      for (auto candidate = children.begin(); candidate != children.end();
           ++candidate) {
        if (retained != children.begin() &&
            (retained - 1)->canonical.key_words ==
                candidate->canonical.key_words) {
          ++result.duplicate_children;
          continue;
        }
        if (retained != candidate) *retained = std::move(*candidate);
        ++retained;
      }
      children.erase(retained, children.end());
      current = std::move(children);
    }
    return false;
  };

  // The certified q=8 primitive-tensor result says every nondegenerate
  // q=8 tensor has a nondegenerate q=7 hyperplane. Reuse the complete q=7
  // seeds for both the k=8 and k=7 sectors.
  auto dimension_seven_seeds = seed_nodes(7);
  if (!dimension_seven_seeds.empty() &&
      sector_has_witness(dimension_seven_seeds, 8)) {
    result.has_isotropic_subspace = true;
    result.sector_witness_found[8] = true;
    return result;
  }
  if (!dimension_seven_seeds.empty() &&
      sector_has_witness(std::move(dimension_seven_seeds), 7)) {
    result.has_isotropic_subspace = true;
    result.sector_witness_found[7] = true;
    return result;
  }

  for (std::uint32_t intrinsic_dimension = 6;
       intrinsic_dimension >= 1; --intrinsic_dimension) {
    auto seeds = seed_nodes(intrinsic_dimension);
    if (!seeds.empty() &&
        sector_has_witness(std::move(seeds), intrinsic_dimension)) {
      result.has_isotropic_subspace = true;
      result.sector_witness_found[intrinsic_dimension] = true;
      return result;
    }
  }
  auto zero_tensor_seed = seed_nodes(0);
  if (sector_has_witness(std::move(zero_tensor_seed), 0)) {
    result.has_isotropic_subspace = true;
    result.sector_witness_found[0] = true;
  }
  return result;
}

}  // namespace utsp
