#include "utsp/origin_orbits.hpp"

#include "utsp/native_core.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace utsp {
namespace {

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

[[nodiscard]] bool parity(std::uint32_t left, std::uint32_t right) {
  return (std::popcount(left & right) & 1U) != 0;
}

[[nodiscard]] std::vector<std::uint32_t> affine_weight_colors(
    std::span<const std::uint32_t> support,
    std::uint32_t ambient_dimension) {
  const auto dual_dimension = ambient_dimension + 1;
  std::vector<std::uint32_t> colors(std::size_t{1} << dual_dimension, 0);
  const auto constant_bit = std::uint32_t{1} << ambient_dimension;
  for (std::uint32_t function = 1; function < colors.size(); ++function) {
    const auto linear = function & (constant_bit - 1);
    const bool constant = (function & constant_bit) != 0;
    std::uint32_t weight = 0;
    for (const auto point : support) {
      weight += static_cast<std::uint32_t>(parity(linear, point) != constant);
    }
    colors[function] = weight;
  }
  return colors;
}

using RefinedColor = std::array<std::uint64_t, 2>;
using MarkedVectorColor =
    std::tuple<
        std::uint32_t,
        std::uint64_t,
        std::uint64_t,
        std::uint64_t,
        bool,
        std::uint64_t,
        std::uint64_t>;

struct MarkedRefinementTable {
  std::size_t dual_count = 0;
  std::vector<std::uint64_t> values;

  [[nodiscard]] std::span<const std::uint64_t> row(
      std::uint32_t origin) const {
    if (values.empty()) return {};
    return std::span<const std::uint64_t>(
        values.data() + static_cast<std::size_t>(origin) * dual_count,
        dual_count);
  }
};

[[nodiscard]] std::uint64_t relation_value(
    std::span<const std::uint64_t> relation,
    std::uint32_t vector) {
  return relation.empty() ? 0 : relation[vector];
}

[[nodiscard]] RefinedColor second_color_value(
    std::span<const RefinedColor> colors,
    std::uint32_t vector) {
  return colors.empty() ? RefinedColor{} : colors[vector];
}

[[nodiscard]] std::uint64_t modular_power(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint64_t modulus) {
  std::uint64_t result = 1;
  while (exponent != 0) {
    if ((exponent & 1U) != 0) result = result * base % modulus;
    base = base * base % modulus;
    exponent >>= 1U;
  }
  return result;
}

void modular_walsh_transform(
    std::span<std::uint64_t> values,
    std::uint64_t modulus,
    bool inverse) {
  for (std::size_t stride = 1; stride < values.size(); stride <<= 1U) {
    for (std::size_t block = 0; block < values.size(); block += 2 * stride) {
      for (std::size_t offset = 0; offset < stride; ++offset) {
        const auto left = values[block + offset];
        const auto right = values[block + stride + offset];
        const auto sum = left + right;
        values[block + offset] = sum >= modulus ? sum - modulus : sum;
        values[block + stride + offset] =
            left >= right ? left - right : left + modulus - right;
      }
    }
  }
  if (!inverse) return;
  const auto inverse_size = modular_power(values.size(), modulus - 2, modulus);
  for (auto& value : values) value = value * inverse_size % modulus;
}

[[nodiscard]] std::vector<RefinedColor> xor_pair_fingerprints(
    std::span<const std::uint32_t> labels,
    std::uint64_t seed) {
  constexpr std::array<std::uint64_t, 2> moduli{
      1'000'000'007ULL,
      1'000'000'009ULL,
  };
  constexpr std::size_t rank = 4;
  const auto label_count =
      *std::max_element(labels.begin(), labels.end()) + 1;
  std::vector<RefinedColor> result(labels.size());
  for (std::size_t channel = 0; channel < moduli.size(); ++channel) {
    const auto modulus = moduli[channel];
    std::vector<std::uint64_t> left_values(labels.size());
    std::vector<std::uint64_t> right_values(labels.size());
    std::vector<std::uint64_t> left_by_label(label_count);
    std::vector<std::uint64_t> right_by_label(label_count);
    for (std::size_t component = 0; component < rank; ++component) {
      const auto component_seed = splitmix64(
          seed ^ (std::uint64_t{channel} << 48U) ^
          (std::uint64_t{component} << 32U));
      for (std::uint32_t label = 0; label < label_count; ++label) {
        left_by_label[label] =
            splitmix64(component_seed ^ label) % modulus;
        right_by_label[label] =
            splitmix64(component_seed ^ 0xd6e8feb86659fd93ULL ^ label) %
            modulus;
      }
      for (std::size_t index = 0; index < labels.size(); ++index) {
        left_values[index] = left_by_label[labels[index]];
        right_values[index] = right_by_label[labels[index]];
      }
      modular_walsh_transform(left_values, modulus, false);
      modular_walsh_transform(right_values, modulus, false);
      for (std::size_t index = 0; index < labels.size(); ++index) {
        left_values[index] = left_values[index] * right_values[index] % modulus;
      }
      modular_walsh_transform(left_values, modulus, true);
      for (std::size_t difference = 0; difference < labels.size(); ++difference) {
        result[difference][channel] += left_values[difference];
        if (result[difference][channel] >= modulus) {
          result[difference][channel] -= modulus;
        }
      }
    }
  }
  return result;
}

[[nodiscard]] MarkedVectorColor marked_vector_color(
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors,
    std::span<const std::uint64_t> relation,
    std::span<const RefinedColor> second_colors,
    std::uint32_t point,
    std::uint32_t vector) {
  const auto second = second_color_value(second_colors, vector);
  return MarkedVectorColor{
      colors[vector],
      refined_colors[vector][0],
      refined_colors[vector][1],
      relation_value(relation, vector),
      parity(point, vector),
      second[0],
      second[1],
  };
}

[[nodiscard]] std::vector<RefinedColor> second_marked_refinement(
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors,
    std::span<const std::uint64_t> relation,
    std::uint32_t point) {
  std::map<MarkedVectorColor, std::uint32_t> color_ids;
  std::vector<std::uint32_t> ids(colors.size());
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    const auto key = marked_vector_color(
        colors, refined_colors, relation, {}, point, vector);
    color_ids.emplace(key, 0);
  }
  std::uint32_t next_id = 0;
  for (auto& [key, identifier] : color_ids) {
    (void)key;
    identifier = next_id++;
  }
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    ids[vector] = color_ids.at(marked_vector_color(
        colors, refined_colors, relation, {}, point, vector));
  }
  return xor_pair_fingerprints(ids, 0x6a09e667f3bcc909ULL);
}

[[nodiscard]] std::map<MarkedVectorColor, std::uint32_t>
marked_color_histogram(
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors,
    std::span<const std::uint64_t> relation,
    std::span<const RefinedColor> second_colors,
    std::uint32_t point) {
  std::map<MarkedVectorColor, std::uint32_t> result;
  for (std::uint32_t vector = 1; vector < colors.size(); ++vector) {
    ++result[marked_vector_color(
        colors, refined_colors, relation, second_colors, point, vector)];
  }
  return result;
}

[[nodiscard]] std::vector<RefinedColor> refined_color_fingerprints(
    std::span<const std::uint32_t> colors,
    bool enabled) {
  std::vector<RefinedColor> result(colors.size());
  if (!enabled) return result;
  return xor_pair_fingerprints(colors, 0x452821e638d01377ULL);
}

[[nodiscard]] std::vector<std::uint32_t> span_in_coefficient_order(
    std::span<const std::uint32_t> basis);

struct PrimalVectorLabel {
  std::uint32_t difference_count = 0;
  RefinedColor difference_refinement{};
  std::uint32_t support_profile = 0;
  std::uint32_t triple_count = 0;
  bool in_translated_support = false;

  [[nodiscard]] auto tuple() const {
    return std::tie(
        difference_count,
        difference_refinement,
        support_profile,
        triple_count,
        in_translated_support);
  }
  [[nodiscard]] bool operator<(const PrimalVectorLabel& other) const {
    return tuple() < other.tuple();
  }
  [[nodiscard]] bool operator==(const PrimalVectorLabel& other) const {
    return tuple() == other.tuple();
  }
};

struct PrimalOriginContext {
  std::uint32_t dimension = 0;
  std::uint32_t point_count = 0;
  std::vector<bool> support;
  std::vector<std::uint32_t> difference_counts;
  std::vector<RefinedColor> difference_refinements;
  std::vector<std::uint32_t> support_profile_ids;
  std::vector<std::uint32_t> triple_counts;

  [[nodiscard]] PrimalVectorLabel label(
      std::uint32_t vector,
      std::uint32_t origin) const {
    const auto translated = vector ^ origin;
    return PrimalVectorLabel{
        difference_counts[vector],
        difference_refinements[vector],
        support_profile_ids[translated],
        triple_counts[translated],
        support[translated],
    };
  }
};

[[nodiscard]] PrimalOriginContext primal_origin_context(
    std::span<const std::uint32_t> support,
    std::uint32_t dimension) {
  PrimalOriginContext result;
  result.dimension = dimension;
  result.point_count = std::uint32_t{1} << dimension;
  result.support.resize(result.point_count, false);
  for (const auto point : support) result.support[point] = true;

  result.difference_counts.resize(result.point_count, 0);
  for (const auto left : support) {
    for (const auto right : support) {
      ++result.difference_counts[left ^ right];
    }
  }
  result.difference_refinements = xor_pair_fingerprints(
      result.difference_counts, 0xbb67ae8584caa73bULL);

  using DifferenceProfile = std::array<std::uint8_t, 65>;
  std::vector<DifferenceProfile> profiles(result.point_count);
  std::map<DifferenceProfile, std::uint32_t> profile_ids;
  for (std::uint32_t point = 0; point < result.point_count; ++point) {
    auto& profile = profiles[point];
    for (const auto support_point : support) {
      ++profile.at(result.difference_counts[point ^ support_point]);
    }
    profile_ids.emplace(profile, 0);
  }
  std::uint32_t next_profile = 0;
  for (auto& [profile, identifier] : profile_ids) {
    (void)profile;
    identifier = next_profile++;
  }
  result.support_profile_ids.resize(result.point_count);
  for (std::uint32_t point = 0; point < result.point_count; ++point) {
    result.support_profile_ids[point] = profile_ids.at(profiles[point]);
  }

  result.triple_counts.resize(result.point_count, 0);
  for (std::size_t first = 0; first < support.size(); ++first) {
    for (std::size_t second = first + 1; second < support.size(); ++second) {
      for (std::size_t third = second + 1; third < support.size(); ++third) {
        ++result.triple_counts[
            support[first] ^ support[second] ^ support[third]];
      }
    }
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> greedy_primal_basis(
    const PrimalOriginContext& context,
    std::span<const std::uint32_t> support,
    std::uint32_t source_origin) {
  std::vector<std::uint32_t> translated_source;
  translated_source.reserve(support.size());
  std::map<PrimalVectorLabel, std::uint32_t> support_label_frequencies;
  for (const auto point : support) {
    const auto vector = point ^ source_origin;
    if (vector != 0) {
      translated_source.push_back(vector);
      ++support_label_frequencies[context.label(vector, source_origin)];
    }
  }

  std::vector<std::uint32_t> basis;
  std::vector<std::uint32_t> span{0};
  while (basis.size() < context.dimension) {
    std::uint32_t best = 0;
    auto best_count = std::numeric_limits<std::uint32_t>::max();
    for (const auto candidate : translated_source) {
      if (std::find(span.begin(), span.end(), candidate) != span.end()) continue;
      const auto count =
          support_label_frequencies.at(context.label(candidate, source_origin));
      if (best == 0 || count < best_count ||
          (count == best_count && candidate < best)) {
        best = candidate;
        best_count = count;
      }
    }
    if (best == 0) {
      throw std::logic_error("failed to construct a primal support basis");
    }
    basis.push_back(best);
    const auto old_size = span.size();
    span.reserve(old_size * 2);
    for (std::size_t index = 0; index < old_size; ++index) {
      span.push_back(span[index] ^ best);
    }
  }
  return basis;
}

struct PrimalEquivalenceSearch {
  const PrimalOriginContext& context;
  std::span<const std::uint32_t> source_by_coeff;
  std::span<const std::vector<std::uint32_t>> candidates_by_depth;
  std::uint32_t source_origin = 0;
  std::uint32_t target_origin = 0;
  std::uint64_t* total_nodes = nullptr;
  std::uint64_t node_limit = 0;
  bool aborted = false;
  std::vector<std::uint32_t> image_by_coeff;
  std::optional<std::vector<std::uint32_t>> image_basis;

  [[nodiscard]] bool run(std::uint32_t depth) {
    if (depth == context.dimension) {
      std::vector<std::uint32_t> images;
      images.reserve(context.dimension);
      for (std::uint32_t index = 0; index < context.dimension; ++index) {
        images.push_back(image_by_coeff[std::size_t{1} << index]);
      }
      image_basis = std::move(images);
      return true;
    }
    const auto prior_count = std::size_t{1} << depth;
    for (const auto candidate : candidates_by_depth[depth]) {
      if (node_limit != 0 && *total_nodes >= node_limit) {
        aborted = true;
        return false;
      }
      ++*total_nodes;
      if (std::find(
              image_by_coeff.begin(),
              image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count),
              candidate) !=
          image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count)) {
        continue;
      }
      bool compatible = true;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        const auto source = source_by_coeff[prior | prior_count];
        const auto image = image_by_coeff[prior] ^ candidate;
        if (context.label(source, source_origin) !=
            context.label(image, target_origin)) {
          compatible = false;
          break;
        }
      }
      if (!compatible) continue;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        image_by_coeff[prior | prior_count] = image_by_coeff[prior] ^ candidate;
      }
      if (run(depth + 1)) return true;
      if (aborted) return false;
    }
    return false;
  }
};

struct PrimalOriginTransport {
  std::vector<std::uint32_t> source_dual_basis;
  std::vector<std::uint32_t> image_dual_basis;
};

struct BasisCoordinateSolver {
  std::vector<std::uint32_t> pivot_vectors;
  std::vector<std::uint32_t> pivot_coefficients;

  explicit BasisCoordinateSolver(std::span<const std::uint32_t> basis)
      : pivot_vectors(basis.size(), 0), pivot_coefficients(basis.size(), 0) {
    for (std::uint32_t index = 0; index < basis.size(); ++index) {
      auto value = basis[index];
      auto coefficients = std::uint32_t{1} << index;
      while (value != 0) {
        const auto pivot = static_cast<std::uint32_t>(std::bit_width(value) - 1);
        if (pivot_vectors[pivot] == 0) {
          pivot_vectors[pivot] = value;
          pivot_coefficients[pivot] = coefficients;
          break;
        }
        value ^= pivot_vectors[pivot];
        coefficients ^= pivot_coefficients[pivot];
      }
      if (value == 0) throw std::logic_error("primal basis is singular");
    }
  }

  [[nodiscard]] std::uint32_t coordinates(std::uint32_t value) const {
    std::uint32_t result = 0;
    while (value != 0) {
      const auto pivot = static_cast<std::uint32_t>(std::bit_width(value) - 1);
      if (pivot >= pivot_vectors.size() || pivot_vectors[pivot] == 0) {
        throw std::logic_error("vector lies outside a primal basis span");
      }
      value ^= pivot_vectors[pivot];
      result ^= pivot_coefficients[pivot];
    }
    return result;
  }
};

[[nodiscard]] std::uint32_t basis_combination(
    std::uint32_t coefficients,
    std::span<const std::uint32_t> basis) {
  std::uint32_t result = 0;
  while (coefficients != 0) {
    const auto bit = static_cast<std::uint32_t>(std::countr_zero(coefficients));
    result ^= basis[bit];
    coefficients &= coefficients - 1;
  }
  return result;
}

[[nodiscard]] std::optional<PrimalOriginTransport>
exact_primal_origin_transport(
    const PrimalOriginContext& context,
    std::span<const std::uint32_t> support,
    std::span<const std::uint32_t> source_basis,
    std::span<const std::uint32_t> source_by_coeff,
    std::uint32_t source_origin,
    std::uint32_t target_origin,
    std::uint64_t& total_nodes,
    std::uint64_t node_limit,
    bool& aborted) {
  std::vector<std::uint32_t> translated_target;
  translated_target.reserve(support.size());
  for (const auto point : support) {
    const auto vector = point ^ target_origin;
    if (vector != 0) translated_target.push_back(vector);
  }
  std::vector<std::vector<std::uint32_t>> candidates_by_depth(
      context.dimension);
  for (std::uint32_t depth = 0; depth < context.dimension; ++depth) {
    const auto required = context.label(source_basis[depth], source_origin);
    for (const auto candidate : translated_target) {
      if (context.label(candidate, target_origin) == required) {
        candidates_by_depth[depth].push_back(candidate);
      }
    }
    if (candidates_by_depth[depth].empty()) return std::nullopt;
  }

  PrimalEquivalenceSearch search{
      context,
      source_by_coeff,
      candidates_by_depth,
      source_origin,
      target_origin,
      &total_nodes,
      node_limit,
      false,
      std::vector<std::uint32_t>(context.point_count, 0),
      std::nullopt,
  };
  const bool found = search.run(0);
  aborted = search.aborted;
  if (!found) return std::nullopt;

  const auto& target_basis = *search.image_basis;
  const BasisCoordinateSolver source_solver(source_basis);
  const BasisCoordinateSolver target_solver(target_basis);
  const auto linear_image = [&](std::uint32_t point) {
    return basis_combination(source_solver.coordinates(point), target_basis);
  };
  const auto inverse_linear_image = [&](std::uint32_t point) {
    return basis_combination(target_solver.coordinates(point), source_basis);
  };
  const auto translation = linear_image(source_origin) ^ target_origin;
  const auto inverse_translation = inverse_linear_image(translation);
  const auto constant = std::uint32_t{1} << context.dimension;

  PrimalOriginTransport result;
  result.source_dual_basis.reserve(context.dimension + 1);
  result.image_dual_basis.reserve(context.dimension + 1);
  result.source_dual_basis.push_back(constant);
  result.image_dual_basis.push_back(constant);
  for (std::uint32_t coordinate = 0;
       coordinate < context.dimension;
       ++coordinate) {
    result.source_dual_basis.push_back(std::uint32_t{1} << coordinate);
    std::uint32_t image_function = 0;
    for (std::uint32_t target_coordinate = 0;
         target_coordinate < context.dimension;
         ++target_coordinate) {
      if (((inverse_linear_image(std::uint32_t{1} << target_coordinate) >>
            coordinate) &
           1U) != 0) {
        image_function |= std::uint32_t{1} << target_coordinate;
      }
    }
    if (((inverse_translation >> coordinate) & 1U) != 0) {
      image_function |= constant;
    }
    result.image_dual_basis.push_back(image_function);
  }
  return result;
}

[[nodiscard]] MarkedRefinementTable marked_refinement_table(
    std::span<const std::uint32_t> colors,
    std::uint32_t ambient_dimension,
    bool enabled) {
  MarkedRefinementTable result;
  result.dual_count = colors.size();
  if (!enabled) return result;
  const auto point_count = std::size_t{1} << ambient_dimension;
  result.values.resize(point_count * colors.size());
  const auto color_count =
      *std::max_element(colors.begin(), colors.end()) + 1;
  std::vector<std::uint64_t> pair_hash(
      static_cast<std::size_t>(color_count) * color_count);
  for (std::uint32_t left = 0; left < color_count; ++left) {
    for (std::uint32_t right = 0; right < color_count; ++right) {
      pair_hash[static_cast<std::size_t>(left) * color_count + right] =
          splitmix64(
              0x3f84d5b5b5470917ULL ^ (std::uint64_t{left} << 32U) ^ right);
    }
  }
  const auto constant = std::uint32_t{1} << ambient_dimension;
  std::vector<std::uint64_t> transform(colors.size());
  for (std::uint32_t difference = 0; difference < colors.size(); ++difference) {
    for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
      transform[vector] = pair_hash[
          static_cast<std::size_t>(colors[vector]) * color_count +
          colors[vector ^ difference]];
    }
    for (std::size_t stride = 1; stride < transform.size(); stride <<= 1U) {
      for (std::size_t block = 0; block < transform.size(); block += 2 * stride) {
        for (std::size_t offset = 0; offset < stride; ++offset) {
          const auto left = transform[block + offset];
          const auto right = transform[block + stride + offset];
          transform[block + offset] = left + right;
          transform[block + stride + offset] = left - right;
        }
      }
    }
    for (std::uint32_t origin = 0; origin < point_count; ++origin) {
      result.values[
          static_cast<std::size_t>(origin) * colors.size() + difference] =
          transform[origin | constant];
    }
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> greedy_marked_basis(
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors,
    std::span<const std::uint64_t> source_relation,
    std::span<const std::uint64_t> target_relation,
    std::span<const RefinedColor> source_second_colors,
    std::span<const RefinedColor> target_second_colors,
    std::uint32_t fixed_vector,
    std::uint32_t dimension,
    std::uint32_t source_point,
    std::uint32_t target_point) {
  std::map<MarkedVectorColor, std::uint32_t> marked_counts;
  for (std::size_t vector = 1; vector < colors.size(); ++vector) {
    const auto marked_color = marked_vector_color(
        colors,
        refined_colors,
        target_relation,
        target_second_colors,
        target_point,
        static_cast<std::uint32_t>(vector));
    ++marked_counts[marked_color];
  }
  std::vector<std::uint32_t> basis{fixed_vector};
  std::vector<std::uint32_t> span{0, fixed_vector};
  while (basis.size() < dimension) {
    std::uint32_t best = 0;
    std::vector<std::uint32_t> best_score;
    for (std::uint32_t candidate = 1; candidate < colors.size(); ++candidate) {
      if (std::find(span.begin(), span.end(), candidate) != span.end()) continue;
      std::vector<std::uint32_t> score;
      score.reserve(span.size());
      for (const auto prior : span) {
        const auto vector = candidate ^ prior;
        const auto marked_color = marked_vector_color(
            colors,
            refined_colors,
            source_relation,
            source_second_colors,
            source_point,
            vector);
        score.push_back(marked_counts[marked_color]);
      }
      std::sort(score.begin(), score.end(), std::greater<>());
      if (best == 0 || score < best_score ||
          (score == best_score && candidate < best)) {
        best = candidate;
        best_score = std::move(score);
      }
    }
    if (best == 0) throw std::logic_error("failed to construct a dual basis");
    basis.push_back(best);
    const auto old_size = span.size();
    span.reserve(old_size * 2);
    for (std::size_t index = 0; index < old_size; ++index) {
      span.push_back(span[index] ^ best);
    }
  }
  return basis;
}

[[nodiscard]] std::vector<std::uint32_t> greedy_unmarked_basis(
    std::span<const std::uint32_t> base_labels,
    std::uint32_t fixed_vector,
    std::uint32_t dimension,
    std::uint32_t point) {
  const auto base_label_count =
      *std::max_element(base_labels.begin(), base_labels.end()) + 1;
  std::vector<std::uint32_t> marked_counts(
      static_cast<std::size_t>(base_label_count) * 2, 0);
  const auto label = [&](std::uint32_t vector) {
    return 2 * base_labels[vector] +
           static_cast<std::uint32_t>(parity(point, vector));
  };
  for (std::uint32_t vector = 1; vector < base_labels.size(); ++vector) {
    ++marked_counts[label(vector)];
  }
  std::vector<std::uint32_t> basis{fixed_vector};
  std::vector<std::uint32_t> span{0, fixed_vector};
  while (basis.size() < dimension) {
    std::uint32_t best = 0;
    auto best_count = std::numeric_limits<std::uint32_t>::max();
    for (std::uint32_t candidate = 1; candidate < base_labels.size(); ++candidate) {
      if (std::find(span.begin(), span.end(), candidate) != span.end()) continue;
      const auto count = marked_counts[label(candidate)];
      if (best == 0 || count < best_count ||
          (count == best_count && candidate < best)) {
        best = candidate;
        best_count = count;
      }
    }
    if (best == 0) throw std::logic_error("failed to construct a dual basis");
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

struct MarkedInvariant {
  bool in_support = false;
  std::uint32_t triple_count = 0;
  std::array<std::uint64_t, 4> color_fourier{};

  [[nodiscard]] auto tuple() const {
    return std::tie(in_support, triple_count, color_fourier);
  }
  [[nodiscard]] bool operator<(const MarkedInvariant& other) const {
    return tuple() < other.tuple();
  }
};

[[nodiscard]] std::vector<MarkedInvariant> marked_invariants(
    std::span<const std::uint32_t> support,
    std::uint32_t ambient_dimension,
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors) {
  const auto point_count = std::size_t{1} << ambient_dimension;
  const auto dual_count = colors.size();
  std::vector<MarkedInvariant> result(point_count);
  for (const auto point : support) result[point].in_support = true;
  for (std::size_t first = 0; first < support.size(); ++first) {
    for (std::size_t second = first + 1; second < support.size(); ++second) {
      for (std::size_t third = second + 1; third < support.size(); ++third) {
        ++result[support[first] ^ support[second] ^ support[third]].triple_count;
      }
    }
  }

  constexpr std::array<std::uint64_t, 4> seeds{
      0x243f6a8885a308d3ULL,
      0x13198a2e03707344ULL,
      0xa4093822299f31d0ULL,
      0x082efa98ec4e6c89ULL,
  };
  const auto constant_bit = std::uint32_t{1} << ambient_dimension;
  for (std::size_t channel = 0; channel < seeds.size(); ++channel) {
    std::vector<std::uint64_t> transform(dual_count);
    for (std::size_t function = 0; function < dual_count; ++function) {
      transform[function] = splitmix64(
          seeds[channel] ^ colors[function] ^ refined_colors[function][0] ^
          std::rotl(refined_colors[function][1], 23));
    }
    for (std::size_t stride = 1; stride < dual_count; stride <<= 1U) {
      for (std::size_t block = 0; block < dual_count; block += 2 * stride) {
        for (std::size_t offset = 0; offset < stride; ++offset) {
          const auto left = transform[block + offset];
          const auto right = transform[block + stride + offset];
          transform[block + offset] = left + right;
          transform[block + stride + offset] = left - right;
        }
      }
    }
    for (std::uint32_t point = 0; point < point_count; ++point) {
      result[point].color_fourier[channel] = transform[point | constant_bit];
    }
  }
  return result;
}

struct EquivalenceSearch {
  std::span<const std::uint32_t> source_labels;
  std::span<const std::uint32_t> target_labels;
  std::span<const std::uint32_t> source_basis;
  std::span<const std::uint32_t> source_by_coeff;
  std::span<const std::vector<std::uint32_t>> vectors_by_label;
  std::uint64_t* total_nodes = nullptr;
  std::uint64_t node_limit = 0;
  bool aborted = false;
  std::vector<std::uint32_t> image_by_coeff;
  std::optional<std::vector<std::uint32_t>> image_basis;

  [[nodiscard]] bool run(std::uint32_t depth) {
    if (depth == source_basis.size()) {
      std::vector<std::uint32_t> images;
      images.reserve(source_basis.size());
      for (std::size_t index = 0; index < source_basis.size(); ++index) {
        images.push_back(image_by_coeff[std::size_t{1} << index]);
      }
      image_basis = std::move(images);
      return true;
    }
    const auto prior_count = std::size_t{1} << depth;
    const auto generator = source_by_coeff[prior_count];
    const auto required_label = source_labels[generator];
    if (required_label >= vectors_by_label.size()) return false;
    for (const auto candidate : vectors_by_label[required_label]) {
      if (node_limit != 0 && *total_nodes >= node_limit) {
        aborted = true;
        return false;
      }
      ++*total_nodes;
      if (std::find(
              image_by_coeff.begin(),
              image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count),
              candidate) !=
          image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count)) {
        continue;
      }
      bool compatible = true;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        const auto image = image_by_coeff[prior] ^ candidate;
        const auto source = source_by_coeff[prior | prior_count];
        if (target_labels[image] != source_labels[source]) {
          compatible = false;
          break;
        }
      }
      if (!compatible) continue;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        image_by_coeff[prior | prior_count] = image_by_coeff[prior] ^ candidate;
      }
      if (run(depth + 1)) return true;
      if (aborted) return false;
    }
    return false;
  }
};

struct UnmarkedEquivalenceSearch {
  std::span<const std::uint32_t> base_labels;
  std::span<const std::uint32_t> source_basis;
  std::span<const std::uint32_t> source_by_coeff;
  std::span<const std::vector<std::uint32_t>> vectors_by_base_label;
  std::uint32_t source_point = 0;
  std::uint32_t target_point = 0;
  std::uint64_t* total_nodes = nullptr;
  std::uint64_t node_limit = 0;
  bool aborted = false;
  std::vector<std::uint32_t> image_by_coeff;
  std::optional<std::vector<std::uint32_t>> image_basis;

  [[nodiscard]] bool run(std::uint32_t depth) {
    if (depth == source_basis.size()) {
      std::vector<std::uint32_t> images;
      images.reserve(source_basis.size());
      for (std::size_t index = 0; index < source_basis.size(); ++index) {
        images.push_back(image_by_coeff[std::size_t{1} << index]);
      }
      image_basis = std::move(images);
      return true;
    }
    const auto prior_count = std::size_t{1} << depth;
    const auto generator = source_by_coeff[prior_count];
    const auto required_base = base_labels[generator];
    const auto required_evaluation = parity(source_point, generator);
    for (const auto candidate : vectors_by_base_label[required_base]) {
      if (parity(target_point, candidate) != required_evaluation) continue;
      if (node_limit != 0 && *total_nodes >= node_limit) {
        aborted = true;
        return false;
      }
      ++*total_nodes;
      if (std::find(
              image_by_coeff.begin(),
              image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count),
              candidate) !=
          image_by_coeff.begin() + static_cast<std::ptrdiff_t>(prior_count)) {
        continue;
      }
      bool compatible = true;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        const auto image = image_by_coeff[prior] ^ candidate;
        const auto source = source_by_coeff[prior | prior_count];
        if (base_labels[image] != base_labels[source] ||
            parity(target_point, image) != parity(source_point, source)) {
          compatible = false;
          break;
        }
      }
      if (!compatible) continue;
      for (std::size_t prior = 0; prior < prior_count; ++prior) {
        image_by_coeff[prior | prior_count] = image_by_coeff[prior] ^ candidate;
      }
      if (run(depth + 1)) return true;
      if (aborted) return false;
    }
    return false;
  }
};

[[nodiscard]] std::optional<std::vector<std::uint32_t>> exact_origin_transport(
    std::span<const std::uint32_t> colors,
    std::span<const RefinedColor> refined_colors,
    std::span<const std::uint32_t> unmarked_base_labels,
    std::span<const std::vector<std::uint32_t>> vectors_by_base_label,
    std::span<const std::uint64_t> source_relation,
    std::span<const std::uint64_t> target_relation,
    std::span<const RefinedColor> source_second_colors,
    std::span<const RefinedColor> target_second_colors,
    std::span<const std::uint32_t> source_basis,
    std::span<const std::uint32_t> source_by_coeff,
    std::uint32_t source_point,
    std::uint32_t target_point,
    std::uint64_t& total_nodes,
    std::uint64_t node_limit,
    bool& aborted) {
  if (source_relation.empty() && target_relation.empty() &&
      source_second_colors.empty() && target_second_colors.empty()) {
    UnmarkedEquivalenceSearch search{
        unmarked_base_labels,
        source_basis,
        source_by_coeff,
        vectors_by_base_label,
        source_point,
        target_point,
        &total_nodes,
        node_limit,
        false,
        std::vector<std::uint32_t>(source_by_coeff.size(), 0),
        std::nullopt,
    };
    search.image_by_coeff[1] = source_basis.front();
    const bool found = search.run(1);
    aborted = search.aborted;
    if (!found) return std::nullopt;
    return search.image_basis;
  }

  std::vector<MarkedVectorColor> source_keys(colors.size());
  std::vector<MarkedVectorColor> target_keys(colors.size());
  std::map<MarkedVectorColor, std::uint32_t> identifiers;
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    source_keys[vector] = marked_vector_color(
        colors,
        refined_colors,
        source_relation,
        source_second_colors,
        source_point,
        vector);
    target_keys[vector] = marked_vector_color(
        colors,
        refined_colors,
        target_relation,
        target_second_colors,
        target_point,
        vector);
    identifiers.emplace(source_keys[vector], 0);
    identifiers.emplace(target_keys[vector], 0);
  }
  std::uint32_t next_identifier = 0;
  for (auto& [key, identifier] : identifiers) {
    (void)key;
    identifier = next_identifier++;
  }
  std::vector<std::uint32_t> source_labels(colors.size());
  std::vector<std::uint32_t> target_labels(colors.size());
  std::vector<std::vector<std::uint32_t>> by_label(identifiers.size());
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    source_labels[vector] = identifiers.at(source_keys[vector]);
    target_labels[vector] = identifiers.at(target_keys[vector]);
    if (vector != 0) by_label[target_labels[vector]].push_back(vector);
  }
  EquivalenceSearch search{
      source_labels,
      target_labels,
      source_basis,
      source_by_coeff,
      by_label,
      &total_nodes,
      node_limit,
      false,
      std::vector<std::uint32_t>(source_by_coeff.size(), 0),
      std::nullopt,
  };
  search.image_by_coeff[1] = source_basis.front();
  const bool found = search.run(1);
  aborted = search.aborted;
  if (!found) return std::nullopt;
  return search.image_basis;
}

}  // namespace

static AffineOriginOrbitResult affine_origin_orbits_impl(
    std::span<const std::uint32_t> support_points,
    std::uint32_t ambient_dimension,
    std::uint64_t node_limit,
    bool use_refined_colors,
    bool use_marked_relations,
    bool use_primal_search) {
  if (ambient_dimension < 1 || ambient_dimension > 16) {
    throw std::invalid_argument("origin orbit search currently supports m=1,...,16");
  }
  if (support_points.empty() || support_points.size() > 64) {
    throw std::invalid_argument("origin orbit support must have 1,...,64 points");
  }
  const auto point_count = std::uint32_t{1} << ambient_dimension;
  std::unordered_set<std::uint32_t> distinct;
  std::vector<Mask> affine_differences;
  affine_differences.reserve(support_points.size() - 1);
  for (const auto point : support_points) {
    if (point >= point_count || !distinct.insert(point).second) {
      throw std::invalid_argument("origin orbit support is invalid");
    }
    if (point != support_points.front()) {
      affine_differences.push_back(point ^ support_points.front());
    }
  }
  if (gf2_rank(affine_differences) != ambient_dimension) {
    throw std::invalid_argument("origin orbit support does not affinely span");
  }

  AffineOriginOrbitResult result;
  result.ambient_dimension = ambient_dimension;
  result.used_refined_colors = use_refined_colors;
  result.used_marked_relations = use_marked_relations;
  std::optional<PrimalOriginContext> primal_context;
  if (use_primal_search) {
    if (!use_refined_colors || use_marked_relations || ambient_dimension > 16) {
      throw std::invalid_argument("primal origin search has invalid tier options");
    }
    primal_context = primal_origin_context(support_points, ambient_dimension);
    result.used_primal_search = true;
  }
  const auto colors = affine_weight_colors(support_points, ambient_dimension);
  const auto refined_colors =
      refined_color_fingerprints(colors, use_refined_colors);
  using UnmarkedBaseKey =
      std::tuple<std::uint32_t, std::uint64_t, std::uint64_t>;
  std::map<UnmarkedBaseKey, std::uint32_t> unmarked_identifiers;
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    unmarked_identifiers.emplace(
        UnmarkedBaseKey{
            colors[vector],
            refined_colors[vector][0],
            refined_colors[vector][1]},
        0);
  }
  std::uint32_t next_unmarked_identifier = 0;
  for (auto& [key, identifier] : unmarked_identifiers) {
    (void)key;
    identifier = next_unmarked_identifier++;
  }
  std::vector<std::uint32_t> unmarked_base_labels(colors.size());
  std::vector<std::vector<std::uint32_t>> vectors_by_base_label(
      unmarked_identifiers.size());
  for (std::uint32_t vector = 0; vector < colors.size(); ++vector) {
    const auto identifier = unmarked_identifiers.at(UnmarkedBaseKey{
        colors[vector],
        refined_colors[vector][0],
        refined_colors[vector][1]});
    unmarked_base_labels[vector] = identifier;
    if (vector != 0) vectors_by_base_label[identifier].push_back(vector);
  }
  const auto marked_refinement = marked_refinement_table(
      colors, ambient_dimension, use_marked_relations);
  auto invariants =
      marked_invariants(
          support_points, ambient_dimension, colors, refined_colors);
  if (use_marked_relations) {
    constexpr std::array<std::uint64_t, 4> seeds{
        0x9216d5d98979fb1bULL,
        0xd1310ba698dfb5acULL,
        0x2ffd72dbd01adfb7ULL,
        0xb8e1afed6a267e96ULL,
    };
    for (std::uint32_t origin = 0; origin < point_count; ++origin) {
      const auto relation = marked_refinement.row(origin);
      for (std::size_t channel = 0; channel < seeds.size(); ++channel) {
        std::uint64_t aggregate = 0;
        for (const auto value : relation) {
          aggregate += splitmix64(seeds[channel] ^ value);
        }
        invariants[origin].color_fourier[channel] ^= aggregate;
      }
    }
  }
  std::map<MarkedInvariant, std::vector<std::uint32_t>> buckets;
  for (std::uint32_t origin = 0; origin < point_count; ++origin) {
    buckets[invariants[origin]].push_back(origin);
  }
  result.invariant_bucket_count = buckets.size();

  const auto constant = std::uint32_t{1} << ambient_dimension;
  std::vector<bool> assigned(point_count, false);

  for (const auto& [invariant, candidates] : buckets) {
    (void)invariant;
    for (const auto representative : candidates) {
      if (assigned[representative]) continue;
      AffineOriginOrbit orbit;
      orbit.representative = representative;
      const auto representative_point = representative | constant;
      const auto representative_relation =
          marked_refinement.row(representative);
      std::vector<std::uint32_t> identity_basis;
      if (primal_context) {
        identity_basis.push_back(constant);
        for (std::uint32_t bit = 0; bit < ambient_dimension; ++bit) {
          identity_basis.push_back(std::uint32_t{1} << bit);
        }
      } else if (representative_relation.empty()) {
        identity_basis = greedy_unmarked_basis(
            unmarked_base_labels,
            constant,
            ambient_dimension + 1,
            representative_point);
      } else {
        identity_basis = greedy_marked_basis(
            colors,
            refined_colors,
            representative_relation,
            representative_relation,
            {},
            {},
            constant,
            ambient_dimension + 1,
            representative_point,
            representative_point);
      }
      const auto source_by_coeff = primal_context
                                       ? std::vector<std::uint32_t>{}
                                       : span_in_coefficient_order(identity_basis);
      std::vector<std::uint32_t> primal_source_basis;
      std::vector<std::uint32_t> primal_source_by_coeff;
      if (primal_context) {
        primal_source_basis = greedy_primal_basis(
            *primal_context, support_points, representative);
        primal_source_by_coeff =
            span_in_coefficient_order(primal_source_basis);
      }
      std::optional<std::vector<RefinedColor>> representative_second;
      std::vector<std::uint32_t> second_identity_basis;
      std::vector<std::uint32_t> second_source_by_coeff;
      orbit.members.push_back(
          OriginTransporter{representative, identity_basis, identity_basis});
      assigned[representative] = true;
      for (const auto candidate : candidates) {
        if (assigned[candidate]) continue;
        ++result.exact_equivalence_tests;
        bool aborted = false;
        if (primal_context) {
          const auto nodes_before = result.exact_search_nodes;
          auto transport = exact_primal_origin_transport(
              *primal_context,
              support_points,
              primal_source_basis,
              primal_source_by_coeff,
              representative,
              candidate,
              result.exact_search_nodes,
              node_limit,
              aborted);
          result.primal_search_nodes +=
              result.exact_search_nodes - nodes_before;
          if (aborted) {
            result.complete = false;
            return result;
          }
          if (!transport) continue;
          assigned[candidate] = true;
          orbit.members.push_back(OriginTransporter{
              candidate,
              std::move(transport->source_dual_basis),
              std::move(transport->image_dual_basis),
          });
          continue;
        }
        bool used_second_order = false;
        constexpr std::uint64_t second_order_node_gate = 10'000;
        auto trial_limit = node_limit;
        bool second_order_trial = false;
        if (use_refined_colors &&
            result.exact_search_nodes <=
                std::numeric_limits<std::uint64_t>::max() -
                    second_order_node_gate) {
          const auto soft_limit =
              result.exact_search_nodes + second_order_node_gate;
          if (node_limit == 0 || soft_limit < node_limit) {
            trial_limit = soft_limit;
            second_order_trial = true;
          }
        }
        auto images = exact_origin_transport(
            colors,
            refined_colors,
            unmarked_base_labels,
            vectors_by_base_label,
            representative_relation,
            marked_refinement.row(candidate),
            {},
            {},
            identity_basis,
            source_by_coeff,
            representative_point,
            candidate | constant,
            result.exact_search_nodes,
            trial_limit,
            aborted);
        if (aborted) {
          if (!second_order_trial ||
              (node_limit != 0 && result.exact_search_nodes >= node_limit)) {
            result.complete = false;
            return result;
          }
          ++result.second_order_refinements;
          used_second_order = true;
          if (!representative_second) {
            representative_second = second_marked_refinement(
                colors,
                refined_colors,
                representative_relation,
                representative_point);
            second_identity_basis = greedy_marked_basis(
                colors,
                refined_colors,
                representative_relation,
                representative_relation,
                *representative_second,
                *representative_second,
                constant,
                ambient_dimension + 1,
                representative_point,
                representative_point);
            second_source_by_coeff =
                span_in_coefficient_order(second_identity_basis);
          }
          const auto candidate_relation = marked_refinement.row(candidate);
          const auto candidate_point = candidate | constant;
          const auto candidate_second = second_marked_refinement(
              colors,
              refined_colors,
              candidate_relation,
              candidate_point);
          if (marked_color_histogram(
                  colors,
                  refined_colors,
                  representative_relation,
                  *representative_second,
                  representative_point) !=
              marked_color_histogram(
                  colors,
                  refined_colors,
                  candidate_relation,
                  candidate_second,
                  candidate_point)) {
            continue;
          }
          aborted = false;
          images = exact_origin_transport(
              colors,
              refined_colors,
              unmarked_base_labels,
              vectors_by_base_label,
              representative_relation,
              candidate_relation,
              *representative_second,
              candidate_second,
              second_identity_basis,
              second_source_by_coeff,
              representative_point,
              candidate_point,
              result.exact_search_nodes,
              node_limit,
              aborted);
          if (aborted) {
            result.complete = false;
            return result;
          }
        }
        if (!images) continue;
        assigned[candidate] = true;
        orbit.members.push_back(OriginTransporter{
            candidate,
            used_second_order ? second_identity_basis : identity_basis,
            std::move(*images),
        });
      }
      result.orbits.push_back(std::move(orbit));
    }
  }
  if (std::find(assigned.begin(), assigned.end(), false) != assigned.end()) {
    throw std::logic_error("origin orbits did not cover the affine space");
  }
  return result;
}

AffineOriginOrbitResult affine_origin_orbits(
    std::span<const std::uint32_t> support_points,
    std::uint32_t ambient_dimension,
    std::uint64_t node_limit) {
  if (ambient_dimension >= 11) {
    constexpr std::uint64_t primal_node_gate = 10'000'000;
    const auto primal_limit = node_limit == 0
                                  ? primal_node_gate
                                  : std::min(node_limit, primal_node_gate);
    auto primal = affine_origin_orbits_impl(
        support_points,
        ambient_dimension,
        primal_limit,
        true,
        false,
        true);
    const auto primal_nodes = primal.exact_search_nodes;
    primal.primal_search_nodes = primal_nodes;
    primal.unmarked_refinement_search_nodes = primal_nodes;
    if (primal.complete ||
        (node_limit != 0 && node_limit <= primal_node_gate)) {
      return primal;
    }

    const auto refined_limit = node_limit == 0
                                   ? std::uint64_t{0}
                                   : node_limit - primal_nodes;
    auto refined = affine_origin_orbits_impl(
        support_points,
        ambient_dimension,
        refined_limit,
        true,
        false,
        false);
    const auto refined_nodes = refined.exact_search_nodes;
    refined.refined_dual_search_nodes = refined_nodes;
    refined.primal_search_nodes = primal_nodes;
    refined.unmarked_refinement_search_nodes = primal_nodes + refined_nodes;
    refined.exact_search_nodes += primal_nodes;
    refined.exact_equivalence_tests += primal.exact_equivalence_tests;
    refined.second_order_refinements += primal.second_order_refinements;
    refined.used_primal_search = true;
    return refined;
  }

  if (ambient_dimension <= 10) {
    constexpr std::uint64_t refined_dual_node_gate = 100'000;
    const auto refined_limit = node_limit == 0
                                   ? refined_dual_node_gate
                                   : std::min(node_limit, refined_dual_node_gate);
    auto refined = affine_origin_orbits_impl(
        support_points,
        ambient_dimension,
        refined_limit,
        true,
        false,
        false);
    const auto refined_nodes = refined.exact_search_nodes;
    refined.unmarked_refinement_search_nodes = refined_nodes;
    refined.refined_dual_search_nodes = refined_nodes;
    if (refined.complete ||
        (node_limit != 0 && node_limit <= refined_dual_node_gate)) {
      return refined;
    }

    const auto after_refined = node_limit == 0
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : node_limit - refined_nodes;
    constexpr std::uint64_t primal_node_gate = 10'000'000;
    const auto primal_limit = std::min(after_refined, primal_node_gate);
    auto primal = affine_origin_orbits_impl(
        support_points,
        ambient_dimension,
        primal_limit,
        true,
        false,
        true);
    const auto primal_nodes = primal.exact_search_nodes;
    primal.refined_dual_search_nodes = refined_nodes;
    primal.primal_search_nodes = primal_nodes;
    primal.unmarked_refinement_search_nodes = refined_nodes + primal_nodes;
    primal.exact_search_nodes += refined_nodes;
    primal.exact_equivalence_tests += refined.exact_equivalence_tests;
    primal.second_order_refinements += refined.second_order_refinements;
    if (primal.complete ||
        (node_limit != 0 && primal.exact_search_nodes >= node_limit)) {
      return primal;
    }

    const auto marked_limit = node_limit == 0
                                  ? std::uint64_t{0}
                                  : node_limit - primal.exact_search_nodes;
    auto marked = affine_origin_orbits_impl(
        support_points,
        ambient_dimension,
        marked_limit,
        true,
        true,
        false);
    marked.refined_dual_search_nodes = refined_nodes;
    marked.primal_search_nodes = primal_nodes;
    marked.unmarked_refinement_search_nodes = refined_nodes + primal_nodes;
    marked.exact_search_nodes += primal.exact_search_nodes;
    marked.exact_equivalence_tests += primal.exact_equivalence_tests;
    marked.second_order_refinements += primal.second_order_refinements;
    marked.used_primal_search = true;
    return marked;
  }

  constexpr std::uint64_t coarse_node_gate = 1'000'000;
  const auto coarse_limit =
      node_limit == 0 ? coarse_node_gate : std::min(node_limit, coarse_node_gate);
  auto coarse = affine_origin_orbits_impl(
      support_points, ambient_dimension, coarse_limit, false, false, false);
  if (coarse.complete || (node_limit != 0 && node_limit <= coarse_node_gate)) {
    return coarse;
  }

  const auto after_coarse_limit = node_limit == 0
                                      ? std::numeric_limits<std::uint64_t>::max()
                                      : node_limit - coarse.exact_search_nodes;
  const std::uint64_t unmarked_refinement_node_gate =
      ambient_dimension <= 10 ? 10'000'000 : 50'000'000;
  const auto unmarked_limit =
      std::min(after_coarse_limit, unmarked_refinement_node_gate);
  auto unmarked = affine_origin_orbits_impl(
      support_points, ambient_dimension, unmarked_limit, true, false, false);
  const auto unmarked_nodes = unmarked.exact_search_nodes;
  unmarked.coarse_search_nodes = coarse.exact_search_nodes;
  unmarked.unmarked_refinement_search_nodes = unmarked_nodes;
  unmarked.exact_search_nodes += coarse.exact_search_nodes;
  unmarked.primal_search_nodes += coarse.primal_search_nodes;
  unmarked.exact_equivalence_tests += coarse.exact_equivalence_tests;
  unmarked.second_order_refinements += coarse.second_order_refinements;
  if (unmarked.complete ||
      (node_limit != 0 && unmarked.exact_search_nodes >= node_limit)) {
    return unmarked;
  }

  const auto marked_limit = node_limit == 0
                                ? std::uint64_t{0}
                                : node_limit - unmarked.exact_search_nodes;
  auto marked = affine_origin_orbits_impl(
      support_points, ambient_dimension, marked_limit, true, true, false);
  marked.coarse_search_nodes = coarse.exact_search_nodes;
  marked.unmarked_refinement_search_nodes = unmarked_nodes;
  marked.exact_search_nodes += unmarked.exact_search_nodes;
  marked.primal_search_nodes += unmarked.primal_search_nodes;
  marked.exact_equivalence_tests += unmarked.exact_equivalence_tests;
  marked.second_order_refinements += unmarked.second_order_refinements;
  marked.used_primal_search = unmarked.used_primal_search;
  return marked;
}

}  // namespace utsp
