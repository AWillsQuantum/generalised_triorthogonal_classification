#include "utsp/puncture_embedding.hpp"

#include "utsp/native_core.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace utsp {
namespace {

struct VectorHash {
  [[nodiscard]] std::size_t operator()(
      const std::vector<std::uint32_t>& values) const noexcept {
    std::size_t hash = 0x9e3779b97f4a7c15ULL;
    for (const auto value : values) {
      hash ^= static_cast<std::size_t>(value) + 0x9e3779b9U +
              (hash << 6U) + (hash >> 2U);
    }
    return hash;
  }
};

struct EmbeddedPatternKey {
  std::vector<std::uint32_t> direction_basis;
  std::vector<std::uint32_t> translated_pattern;

  [[nodiscard]] bool operator==(const EmbeddedPatternKey&) const = default;
};

struct EmbeddedPatternHash {
  [[nodiscard]] std::size_t operator()(
      const EmbeddedPatternKey& key) const noexcept {
    VectorHash hash_vector;
    auto hash = hash_vector(key.direction_basis);
    hash ^= hash_vector(key.translated_pattern) + 0x9e3779b9U +
            (hash << 6U) + (hash >> 2U);
    return hash;
  }
};

void validate_points(
    std::span<const std::uint32_t> points,
    std::uint32_t dimension,
    bool require_nonempty) {
  if (dimension > 20) {
    throw std::invalid_argument("difference profiles currently support dimension<=20");
  }
  if (require_nonempty && points.empty()) {
    throw std::invalid_argument("point set must not be empty");
  }
  const auto bound = std::uint32_t{1} << dimension;
  std::unordered_set<std::uint32_t> distinct;
  for (const auto point : points) {
    if (point >= bound || !distinct.insert(point).second) {
      throw std::invalid_argument("point set is not distinct or exceeds its space");
    }
  }
}

[[nodiscard]] std::uint32_t affine_dimension(
    std::span<const std::uint32_t> points) {
  if (points.empty()) return 0;
  std::vector<Mask> differences;
  differences.reserve(points.size() - 1);
  for (std::size_t index = 1; index < points.size(); ++index) {
    differences.push_back(points[index] ^ points.front());
  }
  return gf2_rank(differences);
}

[[nodiscard]] std::vector<std::uint32_t> greedy_domain_basis(
    std::span<const std::uint32_t> profile,
    std::span<const std::uint32_t> ambient_color_counts,
    std::uint32_t dimension) {
  std::vector<std::uint32_t> basis;
  std::vector<std::uint32_t> span{0};
  for (std::uint32_t depth = 0; depth < dimension; ++depth) {
    std::uint32_t best = 0;
    std::vector<std::uint32_t> best_score;
    for (std::uint32_t candidate = 1; candidate < profile.size(); ++candidate) {
      if (std::find(span.begin(), span.end(), candidate) != span.end()) continue;
      std::vector<std::uint32_t> score;
      score.reserve(span.size());
      for (const auto prior : span) {
        score.push_back(ambient_color_counts[profile[candidate ^ prior]]);
      }
      std::sort(score.begin(), score.end(), std::greater<>());
      if (best == 0 || score < best_score ||
          (score == best_score && candidate < best)) {
        best = candidate;
        best_score = std::move(score);
      }
    }
    if (best == 0) throw std::logic_error("failed to build a template basis");
    basis.push_back(best);
    const auto old_size = span.size();
    span.reserve(old_size * 2);
    for (std::size_t index = 0; index < old_size; ++index) {
      span.push_back(span[index] ^ best);
    }
  }
  return basis;
}

[[nodiscard]] std::vector<std::uint32_t> span_in_coefficient_order(
    std::span<const std::uint32_t> basis) {
  std::vector<std::uint32_t> values(std::size_t{1} << basis.size(), 0);
  for (std::size_t coefficients = 1; coefficients < values.size(); ++coefficients) {
    const auto bit = static_cast<std::size_t>(std::countr_zero(coefficients));
    values[coefficients] = values[coefficients & (coefficients - 1)] ^ basis[bit];
  }
  return values;
}

[[nodiscard]] std::vector<std::uint32_t> canonical_translation(
    std::span<const std::uint32_t> pattern,
    std::span<const std::uint32_t> direction_space) {
  std::vector<std::uint32_t> best;
  for (const auto shift : direction_space) {
    std::vector<std::uint32_t> translated;
    translated.reserve(pattern.size());
    for (const auto point : pattern) translated.push_back(point ^ shift);
    std::sort(translated.begin(), translated.end());
    if (best.empty() || translated < best) best = std::move(translated);
  }
  return best;
}

[[nodiscard]] bool affine_pattern_equal(
    std::span<const std::uint32_t> bucket,
    std::span<const std::uint32_t> pattern) {
  if (bucket.size() != pattern.size()) return false;
  std::vector<std::uint32_t> expected(bucket.begin(), bucket.end());
  std::sort(expected.begin(), expected.end());
  for (const auto target : bucket) {
    for (const auto source : pattern) {
      const auto shift = target ^ source;
      std::vector<std::uint32_t> translated;
      translated.reserve(pattern.size());
      for (const auto point : pattern) translated.push_back(point ^ shift);
      std::sort(translated.begin(), translated.end());
      if (translated == expected) return true;
    }
  }
  return false;
}

[[nodiscard]] std::vector<std::uint32_t> mask_basis_as_u32(
    std::span<const Mask> basis) {
  std::vector<std::uint32_t> result;
  result.reserve(basis.size());
  for (const auto value : basis) result.push_back(static_cast<std::uint32_t>(value));
  return result;
}

}  // namespace

std::vector<std::uint32_t> pair_difference_profile(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension) {
  validate_points(points, ambient_dimension, false);
  std::vector<std::uint32_t> profile(
      std::size_t{1} << ambient_dimension, 0);
  for (std::size_t left = 0; left < points.size(); ++left) {
    for (std::size_t right = left + 1; right < points.size(); ++right) {
      ++profile[points[left] ^ points[right]];
    }
  }
  return profile;
}

HereditaryPunctureSearchResult find_hereditary_punctures(
    std::span<const std::uint32_t> support_points,
    std::uint32_t ambient_dimension,
    std::span<const std::uint32_t> template_points,
    std::uint32_t template_dimension,
    std::uint64_t node_limit) {
  validate_points(support_points, ambient_dimension, true);
  validate_points(template_points, template_dimension, true);
  if (support_points.size() > 64) {
    throw std::invalid_argument("puncture index masks support at most 64 points");
  }
  if (template_dimension > ambient_dimension) {
    throw std::invalid_argument("template dimension exceeds ambient dimension");
  }
  if (affine_dimension(template_points) != template_dimension) {
    throw std::invalid_argument("template points do not affinely span their space");
  }

  HereditaryPunctureSearchResult result;
  if (template_dimension == 0) {
    for (std::size_t index = 0; index < support_points.size(); ++index) {
      result.puncture_index_masks.push_back(std::uint64_t{1} << index);
    }
    result.complete_embeddings = 1;
    result.unique_direction_spaces = 1;
    result.unique_embedded_patterns = 1;
    result.candidate_cosets = support_points.size();
    return result;
  }

  const auto support_profile =
      pair_difference_profile(support_points, ambient_dimension);
  const auto template_profile =
      pair_difference_profile(template_points, template_dimension);
  const auto maximum_color = *std::max_element(
      template_profile.begin(), template_profile.end());
  std::vector<std::uint32_t> ambient_color_counts(maximum_color + 1, 0);
  for (std::size_t direction = 1; direction < support_profile.size(); ++direction) {
    const auto color = support_profile[direction];
    if (color <= maximum_color) ++ambient_color_counts[color];
  }
  std::vector<std::uint32_t> template_color_counts(maximum_color + 1, 0);
  for (std::size_t direction = 1; direction < template_profile.size(); ++direction) {
    ++template_color_counts[template_profile[direction]];
  }
  for (std::size_t color = 0; color < template_color_counts.size(); ++color) {
    if (ambient_color_counts[color] < template_color_counts[color]) {
      result.histogram_rejected = true;
      return result;
    }
  }

  const auto domain_basis = greedy_domain_basis(
      template_profile, ambient_color_counts, template_dimension);
  const auto domain_by_coeff = span_in_coefficient_order(domain_basis);
  std::vector<std::uint32_t> domain_to_coeff(domain_by_coeff.size(), 0);
  for (std::size_t coefficients = 0; coefficients < domain_by_coeff.size();
       ++coefficients) {
    domain_to_coeff[domain_by_coeff[coefficients]] =
        static_cast<std::uint32_t>(coefficients);
  }

  std::vector<std::vector<std::uint32_t>> ambient_by_color(maximum_color + 1);
  for (std::uint32_t direction = 1; direction < support_profile.size(); ++direction) {
    const auto color = support_profile[direction];
    if (color <= maximum_color) ambient_by_color[color].push_back(direction);
  }

  std::vector<std::uint32_t> image_by_coeff(domain_by_coeff.size(), 0);
  std::unordered_set<std::vector<std::uint32_t>, VectorHash> direction_spaces;
  std::unordered_set<EmbeddedPatternKey, EmbeddedPatternHash> embedded_patterns;
  std::unordered_set<std::uint64_t> punctures;
  bool aborted = false;

  std::function<void(std::uint32_t)> extend = [&](std::uint32_t depth) {
    if (aborted) return;
    if (depth == template_dimension) {
      ++result.complete_embeddings;
      std::vector<Mask> image_basis;
      image_basis.reserve(template_dimension);
      for (std::uint32_t index = 0; index < template_dimension; ++index) {
        image_basis.push_back(image_by_coeff[std::size_t{1} << index]);
      }
      const auto reduced_basis = rref_basis(image_basis);
      const auto direction_key = mask_basis_as_u32(reduced_basis);
      direction_spaces.insert(direction_key);

      std::vector<std::uint32_t> image_pattern;
      image_pattern.reserve(template_points.size());
      for (const auto point : template_points) {
        image_pattern.push_back(image_by_coeff[domain_to_coeff[point]]);
      }
      const auto translated_pattern =
          canonical_translation(image_pattern, image_by_coeff);
      EmbeddedPatternKey key{direction_key, translated_pattern};
      if (!embedded_patterns.insert(std::move(key)).second) {
        ++result.duplicate_embedded_patterns;
        return;
      }

      std::unordered_map<std::uint32_t, std::vector<std::size_t>> buckets;
      for (std::size_t index = 0; index < support_points.size(); ++index) {
        const auto representative = static_cast<std::uint32_t>(
            reduce_vector(support_points[index], reduced_basis));
        buckets[representative].push_back(index);
      }
      for (const auto& [representative, indices] : buckets) {
        (void)representative;
        if (indices.size() != template_points.size()) continue;
        ++result.candidate_cosets;
        std::vector<std::uint32_t> bucket_points;
        bucket_points.reserve(indices.size());
        for (const auto index : indices) bucket_points.push_back(support_points[index]);
        if (!affine_pattern_equal(bucket_points, image_pattern)) continue;
        std::uint64_t mask = 0;
        for (const auto index : indices) mask |= std::uint64_t{1} << index;
        punctures.insert(mask);
      }
      return;
    }

    const auto prior_count = std::size_t{1} << depth;
    const auto domain_generator = domain_by_coeff[std::size_t{1} << depth];
    const auto color = template_profile[domain_generator];
    for (const auto candidate : ambient_by_color[color]) {
      if (node_limit != 0 && result.search_nodes >= node_limit) {
        aborted = true;
        result.complete = false;
        return;
      }
      ++result.search_nodes;
      if (std::find(
              image_by_coeff.begin(),
              image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count),
              candidate) !=
          image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count)) {
        ++result.profile_rejections;
        continue;
      }
      bool compatible = true;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        const auto domain = domain_by_coeff[prior | prior_count];
        const auto image = image_by_coeff[prior] ^ candidate;
        if (support_profile[image] != template_profile[domain]) {
          compatible = false;
          break;
        }
      }
      if (!compatible) {
        ++result.profile_rejections;
        continue;
      }
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        image_by_coeff[prior | prior_count] = image_by_coeff[prior] ^ candidate;
      }
      extend(depth + 1);
      if (aborted) return;
    }
  };
  extend(0);

  result.unique_direction_spaces = direction_spaces.size();
  result.unique_embedded_patterns = embedded_patterns.size();
  result.puncture_index_masks.assign(punctures.begin(), punctures.end());
  std::sort(
      result.puncture_index_masks.begin(), result.puncture_index_masks.end());
  return result;
}

}  // namespace utsp
