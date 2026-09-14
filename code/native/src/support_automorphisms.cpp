#include "utsp/support_automorphisms.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace utsp {
namespace {

[[nodiscard]] std::uint32_t pivot_of(Mask value) {
  if (value == 0) throw std::invalid_argument("zero vector has no pivot");
  return 63U - static_cast<std::uint32_t>(std::countl_zero(value));
}

[[nodiscard]] bool dot(Mask left, Mask right) {
  return (std::popcount(left & right) & 1U) != 0;
}

[[nodiscard]] std::uint32_t apply_ambient_map(
    std::uint32_t point,
    std::span<const std::uint32_t> basis_images) {
  std::uint32_t image = 0;
  while (point != 0) {
    const auto index = static_cast<std::size_t>(std::countr_zero(point));
    if (index >= basis_images.size()) {
      throw std::invalid_argument("ambient point lies outside map domain");
    }
    image ^= basis_images[index];
    point &= point - 1;
  }
  return image;
}

[[nodiscard]] bool preserves_support(
    std::span<const std::uint32_t> points,
    std::span<const std::uint32_t> basis_images) {
  std::unordered_set<std::uint32_t> support(points.begin(), points.end());
  for (const auto point : points) {
    if (!support.contains(apply_ambient_map(point, basis_images))) return false;
  }
  return true;
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

[[nodiscard]] Mask permute_row(
    Mask row, std::span<const std::uint32_t> point_images) {
  Mask image = 0;
  while (row != 0) {
    const auto index = static_cast<std::size_t>(std::countr_zero(row));
    image |= Mask{1} << point_images[index];
    row &= row - 1;
  }
  return image;
}

[[nodiscard]] LabelSpaceAutomorphism induce_label_automorphism(
    const LabelSpace& label_space,
    std::vector<std::uint32_t> ambient_basis_images) {
  std::unordered_map<std::uint32_t, std::uint32_t> point_indices;
  point_indices.reserve(label_space.points.size());
  for (std::uint32_t index = 0; index < label_space.points.size(); ++index) {
    point_indices.emplace(label_space.points[index], index);
  }
  std::vector<std::uint32_t> point_images(label_space.points.size());
  for (std::uint32_t index = 0; index < label_space.points.size(); ++index) {
    const auto ambient_image = apply_ambient_map(
        label_space.points[index], ambient_basis_images);
    const auto found = point_indices.find(ambient_image);
    if (found == point_indices.end()) {
      throw std::logic_error("support automorphism moved a point outside support");
    }
    point_images[index] = found->second;
  }

  std::vector<Mask> full_basis = label_space.stabilizer_basis;
  full_basis.insert(
      full_basis.end(),
      label_space.quotient_basis.begin(),
      label_space.quotient_basis.end());
  const auto stabilizer_dimension = label_space.stabilizer_basis.size();
  for (const auto row : label_space.stabilizer_basis) {
    const auto coordinates = coordinates_in_basis(
        permute_row(row, point_images), full_basis);
    if (!coordinates || (*coordinates >> stabilizer_dimension) != 0) {
      throw std::logic_error("ambient map did not preserve the stabilizer subspace");
    }
  }

  LabelSpaceAutomorphism result;
  result.ambient_basis_images = std::move(ambient_basis_images);
  result.quotient_basis_images.reserve(label_space.quotient_dimension());
  for (const auto row : label_space.quotient_basis) {
    const auto coordinates = coordinates_in_basis(
        permute_row(row, point_images), full_basis);
    if (!coordinates) {
      throw std::logic_error("ambient map did not preserve the label space");
    }
    result.quotient_basis_images.push_back(
        *coordinates >> stabilizer_dimension);
  }
  if (gf2_rank(result.quotient_basis_images) !=
      label_space.quotient_dimension()) {
    throw std::logic_error("induced label-space map is singular");
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> ambient_identity(
    std::uint32_t dimension) {
  std::vector<std::uint32_t> result(dimension);
  for (std::uint32_t index = 0; index < dimension; ++index) {
    result[index] = std::uint32_t{1} << index;
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> compose_ambient(
    std::span<const std::uint32_t> left,
    std::span<const std::uint32_t> right) {
  std::vector<std::uint32_t> result;
  result.reserve(right.size());
  for (const auto image : right) {
    result.push_back(apply_ambient_map(image, left));
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> inverse_ambient(
    std::span<const std::uint32_t> map) {
  std::vector<Mask> basis(map.begin(), map.end());
  std::vector<std::uint32_t> result;
  result.reserve(map.size());
  for (std::size_t index = 0; index < map.size(); ++index) {
    const auto coordinates = coordinates_in_basis(Mask{1} << index, basis);
    if (!coordinates || *coordinates >= (Mask{1} << map.size())) {
      throw std::logic_error("ambient automorphism is singular");
    }
    result.push_back(static_cast<std::uint32_t>(*coordinates));
  }
  return result;
}

[[nodiscard]] std::uint64_t ambient_group_order(
    std::span<const LabelSpaceAutomorphism> generators,
    std::uint32_t ambient_dimension) {
  if (generators.empty()) return 1;
  auto current = std::vector<LabelSpaceAutomorphism>(
      generators.begin(), generators.end());
  const auto identity = ambient_identity(ambient_dimension);
  std::uint64_t order = 1;
  const auto point_count = std::size_t{1} << ambient_dimension;
  for (std::uint32_t level = 0; level < ambient_dimension; ++level) {
    const auto base = std::uint32_t{1} << level;
    std::vector<std::int32_t> locations(point_count, -1);
    std::vector<std::uint32_t> orbit{base};
    std::vector<std::vector<std::uint32_t>> transporters{identity};
    locations[base] = 0;
    for (std::size_t next = 0; next < orbit.size(); ++next) {
      for (const auto& generator : current) {
        const auto image = apply_ambient_map(
            orbit[next], generator.ambient_basis_images);
        if (locations[image] >= 0) continue;
        locations[image] = static_cast<std::int32_t>(orbit.size());
        orbit.push_back(image);
        transporters.push_back(compose_ambient(
            generator.ambient_basis_images, transporters[next]));
      }
    }
    if (order > std::numeric_limits<std::uint64_t>::max() / orbit.size()) {
      throw std::overflow_error("support automorphism subgroup order overflow");
    }
    order *= orbit.size();
    std::vector<std::vector<std::uint32_t>> inverse_transporters;
    inverse_transporters.reserve(transporters.size());
    for (const auto& transporter : transporters) {
      inverse_transporters.push_back(inverse_ambient(transporter));
    }
    std::map<std::vector<std::uint32_t>, LabelSpaceAutomorphism> next_generators;
    for (std::size_t index = 0; index < orbit.size(); ++index) {
      for (const auto& generator : current) {
        const auto image = apply_ambient_map(
            orbit[index], generator.ambient_basis_images);
        const auto image_index = static_cast<std::size_t>(locations[image]);
        auto stabilizer_map = compose_ambient(
            inverse_transporters[image_index],
            compose_ambient(
                generator.ambient_basis_images, transporters[index]));
        if (stabilizer_map == identity) continue;
        LabelSpaceAutomorphism retained;
        retained.ambient_basis_images = stabilizer_map;
        next_generators.emplace(
            std::move(stabilizer_map), std::move(retained));
      }
    }
    current.clear();
    current.reserve(next_generators.size());
    for (auto& [key, generator] : next_generators) {
      (void)key;
      current.push_back(std::move(generator));
    }
  }
  return order;
}

[[nodiscard]] std::vector<LabelSpaceAutomorphism> reduce_generators(
    std::span<const LabelSpaceAutomorphism> generators,
    std::uint32_t ambient_dimension) {
  std::vector<LabelSpaceAutomorphism> selected;
  std::uint64_t order = 1;
  for (const auto& generator : generators) {
    auto trial = selected;
    trial.push_back(generator);
    const auto trial_order = ambient_group_order(trial, ambient_dimension);
    if (trial_order == order) continue;
    selected.push_back(generator);
    order = trial_order;
  }
  return selected;
}

[[nodiscard]] LabelSpaceAutomorphism compose_automorphisms(
    const LabelSpaceAutomorphism& left,
    const LabelSpaceAutomorphism& right) {
  LabelSpaceAutomorphism result;
  result.ambient_basis_images = compose_ambient(
      left.ambient_basis_images, right.ambient_basis_images);
  result.quotient_basis_images.reserve(right.quotient_basis_images.size());
  for (const auto image : right.quotient_basis_images) {
    result.quotient_basis_images.push_back(
        linear_combination(image, left.quotient_basis_images));
  }
  return result;
}

[[nodiscard]] LabelSpaceAutomorphism inverse_automorphism(
    const LabelSpaceAutomorphism& automorphism) {
  LabelSpaceAutomorphism result;
  result.ambient_basis_images =
      inverse_ambient(automorphism.ambient_basis_images);
  result.quotient_basis_images.reserve(
      automorphism.quotient_basis_images.size());
  for (std::size_t index = 0;
       index < automorphism.quotient_basis_images.size(); ++index) {
    const auto coordinates = coordinates_in_basis(
        Mask{1} << index, automorphism.quotient_basis_images);
    if (!coordinates) {
      throw std::logic_error("quotient automorphism is singular");
    }
    result.quotient_basis_images.push_back(*coordinates);
  }
  return result;
}

[[nodiscard]] LabelSpaceAutomorphism identity_automorphism(
    std::uint32_t ambient_dimension,
    std::uint32_t quotient_dimension) {
  LabelSpaceAutomorphism result;
  result.ambient_basis_images = ambient_identity(ambient_dimension);
  result.quotient_basis_images.reserve(quotient_dimension);
  for (std::uint32_t index = 0; index < quotient_dimension; ++index) {
    result.quotient_basis_images.push_back(Mask{1} << index);
  }
  return result;
}

[[nodiscard]] SupportAutomorphismSubgroup quotient_vector_stabilizer(
    const SupportAutomorphismSubgroup& subgroup,
    Mask vector) {
  if (subgroup.generators.empty()) return subgroup;
  if (subgroup.quotient_dimension >= 25) {
    throw std::invalid_argument(
        "explicit quotient stabilizers require dimension below 25");
  }
  const auto vector_count =
      std::size_t{1} << subgroup.quotient_dimension;
  std::vector<std::int32_t> locations(vector_count, -1);
  std::vector<Mask> orbit{vector};
  std::vector<LabelSpaceAutomorphism> transporters{
      identity_automorphism(
          subgroup.ambient_dimension, subgroup.quotient_dimension)};
  locations[vector] = 0;
  for (std::size_t next = 0; next < orbit.size(); ++next) {
    for (const auto& generator : subgroup.generators) {
      const auto image = linear_combination(
          orbit[next], generator.quotient_basis_images);
      if (locations[image] >= 0) continue;
      locations[image] = static_cast<std::int32_t>(orbit.size());
      orbit.push_back(image);
      transporters.push_back(
          compose_automorphisms(generator, transporters[next]));
    }
  }
  std::vector<LabelSpaceAutomorphism> inverse_transporters;
  inverse_transporters.reserve(transporters.size());
  for (const auto& transporter : transporters) {
    inverse_transporters.push_back(inverse_automorphism(transporter));
  }
  std::map<std::vector<std::uint32_t>, LabelSpaceAutomorphism> candidates;
  const auto identity = ambient_identity(subgroup.ambient_dimension);
  for (std::size_t index = 0; index < orbit.size(); ++index) {
    for (const auto& generator : subgroup.generators) {
      const auto image = linear_combination(
          orbit[index], generator.quotient_basis_images);
      const auto image_index = static_cast<std::size_t>(locations[image]);
      auto stabilizer = compose_automorphisms(
          inverse_transporters[image_index],
          compose_automorphisms(generator, transporters[index]));
      if (stabilizer.ambient_basis_images == identity) continue;
      candidates.emplace(stabilizer.ambient_basis_images, std::move(stabilizer));
    }
  }
  std::vector<LabelSpaceAutomorphism> candidate_generators;
  candidate_generators.reserve(candidates.size());
  for (auto& [key, generator] : candidates) {
    (void)key;
    candidate_generators.push_back(std::move(generator));
  }
  SupportAutomorphismSubgroup result;
  result.ambient_dimension = subgroup.ambient_dimension;
  result.quotient_dimension = subgroup.quotient_dimension;
  result.elementary_generator_count = candidate_generators.size();
  result.generators = reduce_generators(
      candidate_generators, subgroup.ambient_dimension);
  result.subgroup_order = ambient_group_order(
      result.generators, result.ambient_dimension);
  if (subgroup.subgroup_order % orbit.size() != 0 ||
      result.subgroup_order != subgroup.subgroup_order / orbit.size()) {
    throw std::logic_error("Schreier vector stabilizer has the wrong order");
  }
  for (const auto& generator : result.generators) {
    if (linear_combination(vector, generator.quotient_basis_images) != vector) {
      throw std::logic_error("Schreier generator does not stabilize vector");
    }
  }
  return result;
}

struct ReducedOrbitProblem {
  LabelSpace label_space;
  SupportAutomorphismSubgroup subgroup;
};

[[nodiscard]] std::optional<ReducedOrbitProblem> reduce_after_vector(
    const LabelSpace& label_space,
    const SupportAutomorphismSubgroup& stabilizer,
    Mask representative,
    std::uint32_t remaining_dimension) {
  const auto width = label_space.quotient_dimension();
  std::vector<Mask> constraints;
  constraints.reserve(label_space.bilinear_form_rows.size());
  for (const auto& form : label_space.bilinear_form_rows) {
    constraints.push_back(linear_combination(representative, form));
  }
  const auto common_orthogonal = nullspace_basis(constraints, width);
  std::vector<Mask> span_basis{representative};
  std::vector<Mask> quotient_lifts;
  quotient_lifts.reserve(common_orthogonal.size());
  for (const auto vector : common_orthogonal) {
    const auto residual = reduce_vector(vector, span_basis);
    if (residual == 0) continue;
    quotient_lifts.push_back(residual);
    span_basis.push_back(residual);
    span_basis = rref_basis(span_basis);
  }
  if (quotient_lifts.size() < remaining_dimension) return std::nullopt;

  ReducedOrbitProblem result;
  result.label_space.quotient_basis.reserve(quotient_lifts.size());
  for (const auto lift : quotient_lifts) {
    result.label_space.quotient_basis.push_back(
        linear_combination(lift, label_space.quotient_basis));
  }
  result.label_space.bilinear_form_rows.reserve(
      label_space.bilinear_form_rows.size());
  for (const auto& form : label_space.bilinear_form_rows) {
    std::vector<Mask> induced_rows;
    induced_rows.reserve(quotient_lifts.size());
    for (const auto left : quotient_lifts) {
      const auto contraction = linear_combination(left, form);
      Mask induced_row = 0;
      for (std::size_t right = 0; right < quotient_lifts.size(); ++right) {
        if (dot(contraction, quotient_lifts[right])) {
          induced_row |= Mask{1} << right;
        }
      }
      induced_rows.push_back(induced_row);
    }
    result.label_space.bilinear_form_rows.push_back(std::move(induced_rows));
  }

  result.subgroup.ambient_dimension = stabilizer.ambient_dimension;
  result.subgroup.quotient_dimension = quotient_lifts.size();
  result.subgroup.elementary_generator_count = stabilizer.generators.size();
  result.subgroup.subgroup_order = stabilizer.subgroup_order;
  std::vector<Mask> coordinate_basis{representative};
  coordinate_basis.insert(
      coordinate_basis.end(), quotient_lifts.begin(), quotient_lifts.end());
  result.subgroup.generators.reserve(stabilizer.generators.size());
  for (const auto& generator : stabilizer.generators) {
    LabelSpaceAutomorphism induced;
    induced.ambient_basis_images = generator.ambient_basis_images;
    induced.quotient_basis_images.reserve(quotient_lifts.size());
    for (const auto lift : quotient_lifts) {
      const auto image = linear_combination(
          lift, generator.quotient_basis_images);
      const auto coordinates = coordinates_in_basis(image, coordinate_basis);
      if (!coordinates) {
        throw std::logic_error(
            "vector stabilizer did not preserve the orthogonal quotient");
      }
      induced.quotient_basis_images.push_back(*coordinates >> 1U);
    }
    if (gf2_rank(induced.quotient_basis_images) != quotient_lifts.size()) {
      throw std::logic_error("induced orthogonal-quotient map is singular");
    }
    result.subgroup.generators.push_back(std::move(induced));
  }
  return result;
}

void emit_orbit_cover_candidate(
    std::span<const Mask> quotient_rows,
    std::span<const Mask> physical_basis,
    const IsotropicVisitor& visitor) {
  if (!visitor) return;
  auto canonical_rows = rref_basis(quotient_rows);
  std::vector<Mask> label_rows;
  label_rows.reserve(canonical_rows.size());
  for (const auto row : canonical_rows) {
    label_rows.push_back(linear_combination(row, physical_basis));
  }
  visitor(canonical_rows, label_rows);
}

void enumerate_orbit_cover_recursive(
    const LabelSpace& label_space,
    std::uint32_t remaining_dimension,
    const SupportAutomorphismSubgroup& subgroup,
    std::uint32_t orbit_depth,
    std::vector<Mask>& fixed_rows,
    std::span<const Mask> physical_basis,
    IsotropicOrbitCoverSummary& summary,
    const IsotropicVisitor& visitor) {
  if (remaining_dimension == 0) {
    ++summary.candidate_subspaces;
    emit_orbit_cover_candidate(fixed_rows, physical_basis, visitor);
    return;
  }
  if (orbit_depth == 0 || subgroup.generators.empty() ||
      label_space.quotient_dimension() >= 25) {
    summary.candidate_subspaces += enumerate_totally_isotropic_subspaces(
        label_space,
        remaining_dimension,
        [&](std::span<const Mask> local_rows,
            std::span<const Mask> original_rows) {
          (void)local_rows;
          if (!visitor) return;
          std::vector<Mask> completed = fixed_rows;
          completed.insert(
              completed.end(), original_rows.begin(), original_rows.end());
          emit_orbit_cover_candidate(completed, physical_basis, visitor);
        });
    return;
  }

  const auto orbits = quotient_vector_orbits(subgroup);
  summary.vector_orbits += orbits.size();
  for (const auto& orbit : orbits) {
    const auto original_representative = linear_combination(
        orbit.representative, label_space.quotient_basis);
    fixed_rows.push_back(original_representative);
    if (remaining_dimension == 1) {
      ++summary.candidate_subspaces;
      emit_orbit_cover_candidate(fixed_rows, physical_basis, visitor);
      fixed_rows.pop_back();
      continue;
    }
    const auto stabilizer = quotient_vector_stabilizer(
        subgroup, orbit.representative);
    const auto reduced = reduce_after_vector(
        label_space,
        stabilizer,
        orbit.representative,
        remaining_dimension - 1);
    if (reduced) {
      enumerate_orbit_cover_recursive(
          reduced->label_space,
          remaining_dimension - 1,
          reduced->subgroup,
          orbit_depth - 1,
          fixed_rows,
          physical_basis,
          summary,
          visitor);
    }
    fixed_rows.pop_back();
  }
}

}  // namespace

Mask apply_quotient_automorphism(
    Mask vector, std::span<const Mask> basis_images) {
  return linear_combination(vector, basis_images);
}

std::vector<QuotientVectorOrbit> quotient_vector_orbits(
    const SupportAutomorphismSubgroup& subgroup) {
  if (subgroup.quotient_dimension >= 32) {
    throw std::invalid_argument(
        "explicit quotient-vector orbits require dimension below 32");
  }
  const auto vector_count =
      std::uint64_t{1} << subgroup.quotient_dimension;
  std::vector<bool> seen(static_cast<std::size_t>(vector_count), false);
  std::vector<Mask> queue;
  std::vector<QuotientVectorOrbit> result;
  seen[0] = true;
  for (Mask seed = 1; seed < vector_count; ++seed) {
    if (seen[seed]) continue;
    queue.clear();
    queue.push_back(seed);
    seen[seed] = true;
    std::size_t next = 0;
    Mask representative = seed;
    while (next < queue.size()) {
      const auto vector = queue[next++];
      representative = std::min(representative, vector);
      for (const auto& generator : subgroup.generators) {
        const auto image = apply_quotient_automorphism(
            vector, generator.quotient_basis_images);
        if (seen[image]) continue;
        seen[image] = true;
        queue.push_back(image);
      }
    }
    result.push_back(QuotientVectorOrbit{
        representative, static_cast<std::uint64_t>(queue.size())});
  }
  std::sort(
      result.begin(), result.end(),
      [](const auto& left, const auto& right) {
        return left.representative < right.representative;
      });
  return result;
}

IsotropicOrbitCoverSummary enumerate_isotropic_orbit_cover(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    const SupportAutomorphismSubgroup& subgroup,
    std::uint32_t orbit_depth,
    const IsotropicVisitor& visitor) {
  const auto width = label_space.quotient_dimension();
  if (dimension == 0 || dimension > width) {
    throw std::invalid_argument("logical dimension is outside label quotient");
  }
  if (subgroup.quotient_dimension != width) {
    throw std::invalid_argument("automorphism subgroup has the wrong dimension");
  }
  LabelSpace working;
  working.quotient_basis.reserve(width);
  for (std::uint32_t index = 0; index < width; ++index) {
    working.quotient_basis.push_back(Mask{1} << index);
  }
  working.bilinear_form_rows = label_space.bilinear_form_rows;
  IsotropicOrbitCoverSummary summary;
  std::vector<Mask> fixed_rows;
  fixed_rows.reserve(dimension);
  enumerate_orbit_cover_recursive(
      working,
      dimension,
      subgroup,
      std::min(orbit_depth, dimension),
      fixed_rows,
      label_space.quotient_basis,
      summary,
      visitor);
  return summary;
}

SupportAutomorphismSubgroup elementary_support_automorphisms(
    const LabelSpace& label_space) {
  const auto ambient_dimension = label_space.ambient_dimension;
  if (ambient_dimension == 0 || ambient_dimension >= 32) {
    throw std::invalid_argument("unsupported ambient dimension");
  }

  SupportAutomorphismSubgroup result;
  result.ambient_dimension = ambient_dimension;
  result.quotient_dimension = label_space.quotient_dimension();
  std::set<std::vector<std::uint32_t>> seen;
  auto retain = [&](std::vector<std::uint32_t> basis_images) {
    if (!preserves_support(label_space.points, basis_images) ||
        !seen.insert(basis_images).second) {
      return;
    }
    result.generators.push_back(
        induce_label_automorphism(label_space, std::move(basis_images)));
  };

  for (std::uint32_t left = 0; left < ambient_dimension; ++left) {
    for (std::uint32_t right = left + 1; right < ambient_dimension; ++right) {
      std::vector<std::uint32_t> images(ambient_dimension);
      for (std::uint32_t index = 0; index < ambient_dimension; ++index) {
        images[index] = std::uint32_t{1} << index;
      }
      std::swap(images[left], images[right]);
      retain(std::move(images));
    }
  }
  for (std::uint32_t source = 0; source < ambient_dimension; ++source) {
    for (std::uint32_t target = 0; target < ambient_dimension; ++target) {
      if (source == target) continue;
      std::vector<std::uint32_t> images(ambient_dimension);
      for (std::uint32_t index = 0; index < ambient_dimension; ++index) {
        images[index] = std::uint32_t{1} << index;
      }
      images[source] ^= std::uint32_t{1} << target;
      retain(std::move(images));
    }
  }
  result.elementary_generator_count = result.generators.size();
  result.generators = reduce_generators(
      result.generators, result.ambient_dimension);
  result.subgroup_order = ambient_group_order(
      result.generators, result.ambient_dimension);
  return result;
}

}  // namespace utsp
