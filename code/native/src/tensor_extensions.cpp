#include "utsp/tensor_extensions.hpp"

#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace utsp {
namespace {

struct Monomial {
  std::array<std::uint32_t, 3> indices{};
  std::uint32_t degree = 0;
};

[[nodiscard]] std::vector<Monomial> tensor_monomials(std::uint32_t q) {
  std::vector<Monomial> result;
  result.reserve(cubic_tensor_term_count(q));
  for (std::uint32_t first = 0; first < q; ++first) {
    result.push_back(Monomial{{first, 0, 0}, 1});
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      result.push_back(Monomial{{first, second, 0}, 2});
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      for (std::uint32_t third = second + 1; third < q; ++third) {
        result.push_back(Monomial{{first, second, third}, 3});
      }
    }
  }
  if (result.size() != cubic_tensor_term_count(q)) {
    throw std::logic_error("tensor monomial enumeration has the wrong size");
  }
  return result;
}

[[nodiscard]] bool same_monomial(
    const Monomial& left, const Monomial& right) {
  if (left.degree != right.degree) return false;
  for (std::uint32_t index = 0; index < left.degree; ++index) {
    if (left.indices[index] != right.indices[index]) return false;
  }
  return true;
}

[[nodiscard]] std::uint32_t monomial_position(
    std::span<const Monomial> monomials, const Monomial& selected) {
  for (std::uint32_t index = 0; index < monomials.size(); ++index) {
    if (same_monomial(monomials[index], selected)) return index;
  }
  throw std::logic_error("tensor monomial was not found");
}

[[nodiscard]] std::vector<std::uint32_t> extension_positions(
    std::uint32_t parent_q) {
  const auto child_monomials = tensor_monomials(parent_q + 1);
  std::vector<std::uint32_t> positions;
  for (std::uint32_t position = 0; position < child_monomials.size(); ++position) {
    const auto& monomial = child_monomials[position];
    if (std::find(
            monomial.indices.begin(),
            monomial.indices.begin() + monomial.degree,
            parent_q) != monomial.indices.begin() + monomial.degree) {
      positions.push_back(position);
    }
  }
  return positions;
}

[[nodiscard]] std::uint64_t restrict_standard_hyperplane(
    std::uint32_t child_q, std::uint64_t child_signature) {
  const auto parent_q = child_q - 1;
  const auto parent_monomials = tensor_monomials(parent_q);
  const auto child_monomials = tensor_monomials(child_q);
  std::uint64_t result = 0;
  for (std::uint32_t parent_bit = 0;
       parent_bit < parent_monomials.size(); ++parent_bit) {
    const auto child_bit =
        monomial_position(child_monomials, parent_monomials[parent_bit]);
    result |= ((child_signature >> child_bit) & 1U) << parent_bit;
  }
  return result;
}

[[nodiscard]] std::uint32_t extract_extension_state(
    std::uint32_t parent_q, std::uint64_t child_signature) {
  const auto positions = extension_positions(parent_q);
  std::uint32_t result = 0;
  for (std::uint32_t bit = 0; bit < positions.size(); ++bit) {
    result |= static_cast<std::uint32_t>(
                  (child_signature >> positions[bit]) & 1U)
              << bit;
  }
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

[[nodiscard]] std::uint32_t linear_image(
    std::uint32_t coefficients,
    std::span<const std::uint32_t> basis) {
  std::uint32_t result = 0;
  while (coefficients != 0) {
    const auto least = coefficients & (~coefficients + 1U);
    const auto coordinate =
        static_cast<std::uint32_t>(std::countr_zero(least));
    if (coordinate >= basis.size()) {
      throw std::invalid_argument("basis coefficients are out of range");
    }
    result ^= basis[coordinate];
    coefficients ^= least;
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> identity_basis(
    std::uint32_t dimension) {
  std::vector<std::uint32_t> result(dimension);
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    result[coordinate] = std::uint32_t{1} << coordinate;
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> compose_bases(
    std::span<const std::uint32_t> left,
    std::span<const std::uint32_t> right) {
  if (left.size() != right.size()) {
    throw std::invalid_argument("cannot compose bases of different dimensions");
  }
  std::vector<std::uint32_t> result;
  result.reserve(right.size());
  for (const auto vector : right) {
    result.push_back(linear_image(vector, left));
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> inverse_basis(
    std::span<const std::uint32_t> basis) {
  const auto dimension = static_cast<std::uint32_t>(basis.size());
  if (dimension == 0 || dimension > 7 || rank32(basis) != dimension) {
    throw std::invalid_argument("cannot invert the supplied basis");
  }
  std::vector<std::uint32_t> result(dimension, 0);
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    const auto target = std::uint32_t{1} << coordinate;
    for (std::uint32_t coefficients = 0;
         coefficients < (std::uint32_t{1} << dimension); ++coefficients) {
      if (linear_image(coefficients, basis) == target) {
        result[coordinate] = coefficients;
        break;
      }
    }
    if (result[coordinate] == 0) {
      throw std::logic_error("basis inversion did not find a coordinate");
    }
  }
  return result;
}

[[nodiscard]] std::uint64_t pack_basis(
    std::span<const std::uint32_t> basis) {
  const auto dimension = static_cast<std::uint32_t>(basis.size());
  if (dimension == 0 || dimension > 7 || rank32(basis) != dimension) {
    throw std::invalid_argument("cannot pack the supplied basis");
  }
  std::uint64_t result = 0;
  const auto mask = (std::uint64_t{1} << dimension) - 1;
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    if ((basis[coordinate] & ~mask) != 0) {
      throw std::invalid_argument("basis vector is out of range");
    }
    result |= std::uint64_t{basis[coordinate]}
              << (coordinate * dimension);
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> unpack_basis(
    std::uint64_t packed, std::uint32_t dimension) {
  if (dimension == 0 || dimension > 7) {
    throw std::invalid_argument("packed basis dimension is unsupported");
  }
  const auto mask = (std::uint64_t{1} << dimension) - 1;
  std::vector<std::uint32_t> result(dimension);
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    result[coordinate] = static_cast<std::uint32_t>(
        (packed >> (coordinate * dimension)) & mask);
  }
  if (rank32(result) != dimension) {
    throw std::logic_error("unpacked basis is singular");
  }
  return result;
}

[[nodiscard]] TensorExtensionAffineAction induce_action(
    std::uint32_t parent_q,
    std::uint64_t parent_signature,
    std::span<const std::uint32_t> child_basis) {
  const auto child_q = parent_q + 1;
  if (child_basis.size() != child_q || rank32(child_basis) != child_q) {
    throw std::invalid_argument("extension action basis is not invertible");
  }
  const auto base = embed_tensor_parent_signature(parent_q, parent_signature);
  const auto transformed_base =
      transform_cubic_tensor_signature(child_q, base, child_basis);
  if (restrict_standard_hyperplane(child_q, transformed_base) !=
      parent_signature) {
    throw std::invalid_argument(
        "extension action basis does not stabilize the parent");
  }

  TensorExtensionAffineAction action;
  action.offset = extract_extension_state(parent_q, transformed_base);
  const auto positions = extension_positions(parent_q);
  action.columns.reserve(positions.size());
  for (const auto position : positions) {
    const auto transformed = transform_cubic_tensor_signature(
        child_q, base ^ (std::uint64_t{1} << position), child_basis);
    if (restrict_standard_hyperplane(child_q, transformed) !=
        parent_signature) {
      throw std::logic_error("extension action column changed the parent");
    }
    action.columns.push_back(
        extract_extension_state(parent_q, transformed) ^ action.offset);
  }
  if (rank32(action.columns) != action.columns.size()) {
    throw std::logic_error("extension action has singular linear part");
  }
  return action;
}

class PackedAffineAction {
 public:
  explicit PackedAffineAction(const TensorExtensionAffineAction& action)
      : offset_(action.offset),
        tables_((action.columns.size() + kChunkBits - 1) / kChunkBits) {
    for (std::size_t chunk = 0; chunk < tables_.size(); ++chunk) {
      auto& table = tables_[chunk];
      table.fill(0);
      for (std::uint32_t mask = 1; mask < table.size(); ++mask) {
        const auto least = mask & (~mask + 1U);
        const auto local_bit = static_cast<std::uint32_t>(std::countr_zero(least));
        const auto global_bit =
            static_cast<std::uint32_t>(chunk * kChunkBits) + local_bit;
        table[mask] = table[mask ^ least];
        if (global_bit < action.columns.size()) {
          table[mask] ^= action.columns[global_bit];
        }
      }
    }
  }

  [[nodiscard]] std::uint32_t apply(std::uint32_t value) const {
    auto result = offset_;
    for (std::size_t chunk = 0; chunk < tables_.size(); ++chunk) {
      result ^= tables_[chunk][(value >> (chunk * kChunkBits)) & 0xffU];
    }
    return result;
  }

 private:
  static constexpr std::uint32_t kChunkBits = 8;
  std::uint32_t offset_;
  std::vector<std::array<std::uint32_t, 256>> tables_;
};

}  // namespace

std::uint32_t tensor_extension_bit_count(std::uint32_t parent_q) {
  if (parent_q < 1 || parent_q > 6) {
    throw std::invalid_argument("tensor extension supports parent q=1,...,6");
  }
  return 1 + parent_q + parent_q * (parent_q - 1) / 2;
}

std::uint64_t embed_tensor_parent_signature(
    std::uint32_t parent_q, std::uint64_t parent_signature) {
  const auto parent_monomials = tensor_monomials(parent_q);
  if (parent_signature >= (std::uint64_t{1} << parent_monomials.size())) {
    throw std::invalid_argument("parent tensor signature is out of range");
  }
  const auto child_monomials = tensor_monomials(parent_q + 1);
  std::uint64_t result = 0;
  for (std::uint32_t parent_bit = 0;
       parent_bit < parent_monomials.size(); ++parent_bit) {
    if (((parent_signature >> parent_bit) & 1U) == 0) continue;
    result |= std::uint64_t{1}
              << monomial_position(child_monomials, parent_monomials[parent_bit]);
  }
  return result;
}

std::uint64_t combine_tensor_parent_and_extension(
    std::uint32_t parent_q,
    std::uint64_t parent_signature,
    std::uint32_t extension_state) {
  const auto bit_count = tensor_extension_bit_count(parent_q);
  if (extension_state >= (std::uint32_t{1} << bit_count)) {
    throw std::invalid_argument("tensor extension state is out of range");
  }
  auto result = embed_tensor_parent_signature(parent_q, parent_signature);
  const auto positions = extension_positions(parent_q);
  if (positions.size() != bit_count) {
    throw std::logic_error("tensor extension position count is wrong");
  }
  for (std::uint32_t bit = 0; bit < bit_count; ++bit) {
    result |= std::uint64_t{(extension_state >> bit) & 1U} << positions[bit];
  }
  return result;
}

std::uint64_t restrict_tensor_to_standard_hyperplane(
    std::uint32_t child_q, std::uint64_t child_signature) {
  if (child_q < 2 || child_q > 7) {
    throw std::invalid_argument(
        "standard tensor restriction supports child q=2,...,7");
  }
  return restrict_standard_hyperplane(child_q, child_signature);
}

std::uint32_t tensor_extension_state(
    std::uint32_t parent_q, std::uint64_t child_signature) {
  if (parent_q < 1 || parent_q > 6) {
    throw std::invalid_argument(
        "tensor extension extraction supports parent q=1,...,6");
  }
  return extract_extension_state(parent_q, child_signature);
}

TensorExtensionCensus census_tensor_extensions(
    std::uint32_t parent_q,
    std::uint64_t parent_signature,
    std::span<const std::vector<std::uint32_t>> parent_generators,
    std::uint64_t parent_group_order,
    bool retain_orbit_index_by_state,
    bool retain_transport_to_representative_basis) {
  if (parent_group_order == 0 ||
      parent_group_order >
          (std::numeric_limits<std::uint64_t>::max() >> parent_q)) {
    throw std::invalid_argument("parent automorphism order is invalid");
  }

  TensorExtensionCensus result;
  result.parent_logical_qubits = parent_q;
  result.parent_standard_signature = parent_signature;
  result.extension_bit_count = tensor_extension_bit_count(parent_q);
  result.extension_state_count =
      std::uint32_t{1} << result.extension_bit_count;
  result.parent_automorphism_group_order = parent_group_order;
  result.marked_stabilizer_group_order = parent_group_order << parent_q;

  const auto child_q = parent_q + 1;
  std::vector<std::vector<std::uint32_t>> child_generator_bases;
  for (const auto& parent_basis : parent_generators) {
    if (parent_basis.size() != parent_q || rank32(parent_basis) != parent_q) {
      throw std::invalid_argument("parent automorphism basis is invalid");
    }
    if (transform_cubic_tensor_signature(
            parent_q, parent_signature, parent_basis) != parent_signature) {
      throw std::invalid_argument("parent generator does not preserve its tensor");
    }
    auto child_basis = parent_basis;
    child_basis.push_back(std::uint32_t{1} << parent_q);
    child_generator_bases.push_back(child_basis);
    result.actions.push_back(
        induce_action(parent_q, parent_signature, child_basis));
  }
  for (std::uint32_t coordinate = 0; coordinate < parent_q; ++coordinate) {
    std::vector<std::uint32_t> basis(child_q);
    for (std::uint32_t index = 0; index < child_q; ++index) {
      basis[index] = std::uint32_t{1} << index;
    }
    basis.back() ^= std::uint32_t{1} << coordinate;
    child_generator_bases.push_back(basis);
    result.actions.push_back(induce_action(parent_q, parent_signature, basis));
  }
  if (result.actions.empty()) {
    throw std::logic_error("tensor extension action has no generators");
  }

  std::vector<PackedAffineAction> packed_actions;
  packed_actions.reserve(result.actions.size());
  for (const auto& action : result.actions) packed_actions.emplace_back(action);

  std::vector<std::vector<std::uint32_t>> inverse_generator_bases;
  if (retain_transport_to_representative_basis) {
    inverse_generator_bases.reserve(child_generator_bases.size());
    for (const auto& basis : child_generator_bases) {
      inverse_generator_bases.push_back(inverse_basis(basis));
    }
  }

  constexpr auto unowned = std::numeric_limits<std::uint32_t>::max();
  std::vector<std::uint32_t> owner(result.extension_state_count, unowned);
  std::vector<std::uint64_t> transporters;
  const auto packed_identity = pack_basis(identity_basis(child_q));
  if (retain_transport_to_representative_basis) {
    transporters.resize(result.extension_state_count, 0);
  }
  std::vector<std::uint32_t> queue;
  for (std::uint32_t start = 0; start < result.extension_state_count; ++start) {
    if (owner[start] != unowned) continue;
    const auto orbit_index = static_cast<std::uint32_t>(result.orbits.size());
    queue.clear();
    queue.push_back(start);
    owner[start] = orbit_index;
    if (retain_transport_to_representative_basis) {
      transporters[start] = packed_identity;
    }
    std::size_t head = 0;
    while (head < queue.size()) {
      const auto current = queue[head++];
      for (std::size_t generator = 0;
           generator < packed_actions.size(); ++generator) {
        const auto image = packed_actions[generator].apply(current);
        ++result.generator_transitions;
        if (image >= result.extension_state_count) {
          throw std::logic_error("extension action left its state space");
        }
        if (owner[image] == unowned) {
          owner[image] = orbit_index;
          if (retain_transport_to_representative_basis) {
            const auto current_basis =
                unpack_basis(transporters[current], child_q);
            const auto image_basis = compose_bases(
                inverse_generator_bases[generator], current_basis);
            transporters[image] = pack_basis(image_basis);
          }
          queue.push_back(image);
        }
      }
    }
    const auto orbit_size = static_cast<std::uint32_t>(queue.size());
    if (result.marked_stabilizer_group_order % orbit_size != 0) {
      throw std::logic_error("extension orbit size does not divide group order");
    }
    result.orbits.push_back(TensorExtensionOrbit{start, orbit_size});
  }
  const auto mass = std::accumulate(
      result.orbits.begin(), result.orbits.end(), std::uint64_t{0},
      [](std::uint64_t total, const TensorExtensionOrbit& orbit) {
        return total + orbit.orbit_size;
      });
  if (mass != result.extension_state_count) {
    throw std::logic_error("extension orbit census did not cover every state");
  }
  if (retain_transport_to_representative_basis) {
    const auto stride = std::max<std::uint32_t>(
        1, result.extension_state_count / 4096U);
    for (std::uint32_t state = 0; state < result.extension_state_count;
         state += stride) {
      ++result.transporter_replay_checks;
      const auto orbit_index = owner[state];
      const auto representative = result.orbits[orbit_index].representative;
      const auto basis = unpack_basis(transporters[state], child_q);
      const auto transformed = transform_cubic_tensor_signature(
          child_q,
          combine_tensor_parent_and_extension(
              parent_q, parent_signature, state),
          basis);
      if (transformed != combine_tensor_parent_and_extension(
                             parent_q, parent_signature, representative)) {
        ++result.transporter_replay_mismatches;
      }
    }
    if (result.transporter_replay_mismatches != 0) {
      throw std::logic_error("extension transporter replay failed");
    }
  }
  if (retain_orbit_index_by_state || retain_transport_to_representative_basis) {
    result.orbit_index_by_state = std::move(owner);
  }
  if (retain_transport_to_representative_basis) {
    result.transport_to_representative_basis_by_state =
        std::move(transporters);
  }
  return result;
}

}  // namespace utsp
