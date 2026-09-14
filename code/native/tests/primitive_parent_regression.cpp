#include "utsp/marked_code_canonical.hpp"
#include "utsp/native_core.hpp"
#include "utsp/tensor_canonical.hpp"

#include <array>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  try {
    // Eq. (74) of arXiv:2608.09727v1, expanded through three T parities.
    // The logical basis gives CCZ145 CCZ235, with (n,S,d,A_d)=(48,11,3,240).
    constexpr std::array<utsp::Mask, 11> rows{
        0xdb66dbb6d000ULL, 0x6dbb6ddb6000ULL, 0xb98b98b98b98ULL,
        0x770770770770ULL, 0x7a8d58af0000ULL, 0x6dbdb6b6d000ULL,
        0xb6d6dbdb6000ULL, 0x7a87a87a87a8ULL, 0xaf0af0af0af0ULL,
        0xb6db6db6db6dULL, 0xdb6db6db6db6ULL};
    const auto logical = std::span<const utsp::Mask>(rows).first(5);
    std::vector<std::uint32_t> points;
    for (std::uint32_t j = 0; j < 48; ++j) {
      std::uint32_t point = 0;
      for (std::uint32_t i = 0; i < 6; ++i) {
        point |= ((rows[5+i] >> j) & 1U) << i;
      }
      points.push_back(point);
    }
    if (utsp::gf2_rank(rows) != 11 ||
        utsp::cubic_tensor_word(logical) != 0x500000ULL) {
      throw std::logic_error("the positive witness changed");
    }
    const utsp::DistanceContext distance_context(points, 6);
    const auto distance = distance_context.evaluate(logical);
    if (distance.distance != 3 || distance.error_coefficient != 240) {
      throw std::logic_error("the witness distance changed");
    }
    // The support automorphism action requires the full label quotient.
    const auto space = utsp::build_logical_label_space(points, 6, true);
    std::vector<utsp::Mask> witness_quotient;
    for (const auto row : logical) {
      auto augmented = space.stabilizer_basis;
      augmented.insert(augmented.end(), space.quotient_basis.begin(),
                       space.quotient_basis.end());
      const auto width = augmented.size();
      augmented.push_back(row);
      std::vector<utsp::Mask> constraints;
      for (std::uint32_t column = 0; column < 48; ++column) {
        utsp::Mask constraint = 0;
        for (std::size_t i = 0; i < augmented.size(); ++i) {
          constraint |= ((augmented[i] >> column) & 1U) << i;
        }
        constraints.push_back(constraint);
      }
      const auto relations = utsp::nullspace_basis(constraints, width + 1);
      if (relations.size() != 1 || (relations.front() >> width) != 1) {
        throw std::logic_error("witness is not in the full label space");
      }
      witness_quotient.push_back((relations.front() >> space.stabilizer_basis.size()) &
                                ((utsp::Mask{1} << space.quotient_dimension()) - 1));
    }
    const auto witness_key = utsp::canonicalize_marked_code(
        space, witness_quotient).key_words;
    constexpr std::array<std::uint64_t, 2> targets{0x60000ULL,0x60001ULL};
    const auto reference = utsp::enumerate_target_tensor_subspace_orbits(
        space,5,5,targets,3,1);
    if (std::find(reference.canonical_keys.begin(), reference.canonical_keys.end(),
                  witness_key) == reference.canonical_keys.end()) {
      throw std::logic_error("seeded enumeration missed the literature witness orbit");
    }
    bool okay = true;
    const auto compare = [&](const std::string& mode, const utsp::IsotropicOrbitLevel& result) {
      const bool same = result.canonical_keys == reference.canonical_keys &&
                        result.orbit_sizes == reference.orbit_sizes &&
                        result.weighted_subspace_count == reference.weighted_subspace_count;
      std::cout << "{\"mode\":\"" << mode << "\",\"orbits\":"
                << result.representatives.size() << ",\"matches_seeded_serial\":"
                << (same ? "true" : "false") << "}" << std::endl;
      okay &= same;
    };
    compare("target_seed_three", reference);
    compare("target_seed_three_parallel", utsp::enumerate_target_tensor_subspace_orbits(
        space,5,5,targets,3,2));
    if (argc == 2 && std::string(argv[1]) == "--independent") {
      compare("zero_hyperplane", utsp::enumerate_q5_zero_hyperplane_primitive_subspace_orbits(space));
      compare("target_without_seed", utsp::enumerate_target_tensor_subspace_orbits(
          space,5,5,targets,0,1));
    } else if (argc != 1) {
      throw std::invalid_argument("usage: primitive-parent-regression [--independent]");
    }
    return okay ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
