// Private differential checks keep fixed-support keys out of the public API.
#include "../src/marked_code_canonical.cpp"

#include <iostream>
#include <random>

namespace {

utsp::LabelSpace restrict_space(
    const utsp::LabelSpace& space, const std::vector<utsp::Mask>& basis) {
  auto result = space;
  result.quotient_basis.clear();
  result.bilinear_form_rows.clear();
  for (const auto row : basis) {
    result.quotient_basis.push_back(utsp::linear_combination(row, space.quotient_basis));
  }
  for (const auto& form : space.bilinear_form_rows) {
    std::vector<utsp::Mask> induced;
    for (const auto left : basis) {
      utsp::Mask row = 0;
      const auto contraction = utsp::linear_combination(left, form);
      for (std::size_t right = 0; right < basis.size(); ++right) {
        if (std::popcount(contraction & basis[right]) & 1U) row |= utsp::Mask{1} << right;
      }
      induced.push_back(row);
    }
    result.bilinear_form_rows.push_back(std::move(induced));
  }
  return result;
}

void require(bool condition, const char* message) {
  if (!condition) throw std::logic_error(message);
}

void check_line_parent_orbit(
    const utsp::SmallQuotientGroup& group, std::span<const utsp::Mask> rows) {
  const auto node = group.classify(rows);
  const std::array<utsp::Mask, 3> lines{node.rows[0], node.rows[1], node.rows[0] ^ node.rows[1]};
  std::set<utsp::Mask> accepted;
  for (const auto line : lines) {
    if (group.accepts_line_parent(line, node.rows, node.local_key)) accepted.insert(line);
  }
  require(!accepted.empty(), "canonical incidence selector rejected every parent line");
  std::set<utsp::Mask> orbit{*accepted.begin()};
  for (const auto& element : group.stabilizer_generators(node)) {
    orbit.insert(utsp::linear_combination(*accepted.begin(), element));
  }
  require(accepted == orbit, "accepted lines are not exactly one child-stabilizer orbit");
}

std::size_t check_all_line_parent_orbits(const utsp::SmallQuotientGroup& group,
                                        std::uint32_t width) {
  const auto limit = utsp::Mask{1} << width;
  std::set<std::vector<utsp::Mask>> planes;
  for (utsp::Mask first = 1; first < limit; ++first) {
    for (utsp::Mask second = first + 1; second < limit; ++second) {
      const std::array<utsp::Mask, 2> rows{first, second};
      auto reduced = utsp::rref_basis(rows);
      if (planes.insert(reduced).second) check_line_parent_orbit(group, reduced);
    }
  }
  require(planes.size() == (limit - 1) * (limit - 2) / 6,
          "two-space fixture enumeration does not have the Gaussian-binomial count");
  return planes.size();
}

}  // namespace

int main() {
  try {
    utsp::MarkedCodeCanonicalForm synthetic;
    synthetic.automorphism_group_order = 8;
    synthetic.quotient_automorphism_generators = {{2, 1, 4}};
    const auto kernel_example = utsp::SmallQuotientGroup::build(synthetic, 3);
    require(kernel_example && kernel_example->elements().size() == 2,
            "quotient image closure failed");
    require(kernel_example->classify(std::array<utsp::Mask, 1>{1})
                    .canonical.automorphism_group_order == 4,
            "nonfaithful action lost its kernel");
    require(kernel_example->classify(std::array<utsp::Mask, 1>{3})
                    .canonical.automorphism_group_order == 8,
            "fixed quotient line has the wrong full stabilizer");
    require(!utsp::SmallQuotientGroup::build(synthetic, 3, 1),
            "group size cap was not respected");
    std::size_t incidence_checks = check_all_line_parent_orbits(*kernel_example, 3);
    utsp::MarkedCodeCanonicalForm symmetric;
    symmetric.automorphism_group_order = 480;
    symmetric.quotient_automorphism_generators = {{2, 1, 4, 8, 16}, {2, 4, 8, 16, 1}};
    const auto large_kernel_example = utsp::SmallQuotientGroup::build(symmetric, 5);
    require(large_kernel_example && large_kernel_example->elements().size() == 120,
            "explicit S5 image did not close");
    require(!utsp::SmallQuotientGroup::build(symmetric, 5, 64),
            "restricted image-size cap did not fall back");
    const auto fixed_line = large_kernel_example->classify(std::array<utsp::Mask, 1>{31});
    require(fixed_line.canonical.automorphism_group_order == 480 &&
            fixed_line.stabilizer_indices.size() == 120 && fixed_line.small_stabilizer_mask == 0,
            "large nonfaithful action lost its stabilizer or kernel");
    const auto coordinate_line = large_kernel_example->classify(std::array<utsp::Mask, 1>{1});
    require(coordinate_line.canonical.automorphism_group_order == 96 &&
            large_kernel_example->stabilizer_generators(coordinate_line).size() == 23,
            "large explicit stabilizer has the wrong order");
    const auto large_incidence_checks = check_all_line_parent_orbits(*large_kernel_example, 5);
    incidence_checks += large_incidence_checks;
    utsp::MarkedCodeCanonicalForm cap_fixture;
    cap_fixture.automorphism_group_order = 2880;
    cap_fixture.quotient_automorphism_generators = {
        {2,1,4,8,16,32,64,128}, {2,4,8,16,32,1,64,128}, {1,2,4,8,16,32,128,64}};
    const auto cap_group = utsp::SmallQuotientGroup::build(cap_fixture, 8);
    require(cap_group && cap_group->elements().size() == 1440,
            "S6 times S2 image did not close");
    require(!utsp::SmallQuotientGroup::build(cap_fixture, 8, 1024),
            "restricted image cap did not fall back");
    const auto cap_line = cap_group->classify(std::array<utsp::Mask, 1>{1});
    require(cap_line.canonical.automorphism_group_order == 480 &&
            cap_line.stabilizer_indices.size() == 240,
            "large-cap nonfaithful stabilizer order changed");
    for (utsp::Mask first = 1; first <= 32; ++first) {
      check_line_parent_orbit(*cap_group, std::array<utsp::Mask, 2>{first, 128});
      ++incidence_checks;
    }

    // Fixed supports with independently checked automorphism orders.
    const std::vector<std::vector<std::uint32_t>> fixtures{
      {0,6,63,57,15,17,23,39,46,40,20,34,30,71,113,77,116,92,83,94,98,100,88,
       104,110,82,91,87,107,109,81,79,117,73,69,121,127,67,76,112,64,31,35,41,
       47,19,32,22,16,8,4,2,1},
      {0,15,12,20,27,23,24,117,122,118,121,103,100,107,104,18,29,112,127,115,
       124,110,98,97,109,10,6,83,92,80,95,67,76,79,64,36,41,52,59,57,54,39,42,
       33,45,44,32,31,16,8,4,2,1},
      {204,202,198,3,195,207,5,197,18,20,212,221,34,226,45,237,43,39,235,231,
       225,46,40,232,48,54,58,243,57,249,157,131,158,128,93,67,94,64,49,55,59,
       50,47,32,31,22,21,26,25,16,8,4,2,1},
      {48,235,219,46,30,245,197,192,240,53,198,246,3,216,232,237,221,54,226,
       210,57,18,34,252,204,201,249,39,23,12,63,228,212,10,174,158,176,128,110,
       94,112,64,62,61,19,35,32,16,11,8,7,4,2,1},
      {0,80,75,85,126,125,101,102,116,114,105,111,78,117,118,107,104,122,98,
       100,124,84,87,81,82,69,93,88,64,56,59,38,37,61,62,35,32,19,11,14,22,
       25,26,28,31,13,21,16,8,7,4,2,1}};
    constexpr std::array<std::uint64_t, 5> orders{1, 32, 8, 64, 384};
    std::mt19937_64 random(2026090601);
    std::size_t comparisons = 0;
    std::size_t enumerations = 0;
    std::size_t early_restriction_comparisons = 0;
    for (std::size_t fixture = 0; fixture < fixtures.size(); ++fixture) {
      const auto space = utsp::build_logical_label_space(
          fixtures[fixture], fixture < 2 || fixture == 4 ? 7 : 8, true);
      const auto root = utsp::canonicalize_marked_code(space, {});
      require(root.automorphism_group_order == orders[fixture],
              "support automorphism order differs from the fixture");
      const auto group = utsp::SmallQuotientGroup::build(root, space.quotient_dimension());
      require(group.has_value(), "small fixture group did not close");
      for (std::uint32_t dimension = 1; dimension <= 4; ++dimension) {
        for (std::size_t sample = 0; sample < 4; ++sample) {
          std::vector<utsp::Mask> rows;
          while (rows.size() < dimension) {
            rows.push_back(random() & ((utsp::Mask{1} << space.quotient_dimension()) - 1));
            rows = utsp::rref_basis(rows);
          }
          const auto local = group->classify(rows);
          const auto canonical = utsp::canonicalize_marked_code(space, rows);
          require(local.canonical.key_words.empty(), "local key leaked into global key");
          require(local.canonical.quotient_automorphism_generators.empty(),
                  "small-group node retains redundant action matrices");
          std::vector<std::vector<utsp::Mask>> expected_generators;
          utsp::Mask expected_membership = 0;
          std::vector<std::uint32_t> expected_indices;
          for (std::size_t index = 0; index < group->elements().size(); ++index) {
            std::vector<utsp::Mask> image;
            for (const auto row : local.rows) {
              image.push_back(utsp::linear_combination(row, group->elements()[index]));
            }
            if (utsp::rref_basis(image) != local.rows) continue;
            if (group->elements().size() <= 64) expected_membership |= utsp::Mask{1} << index;
            else expected_indices.push_back(static_cast<std::uint32_t>(index));
            if (index != 0) expected_generators.push_back(group->elements()[index]);
          }
          require(local.small_stabilizer_mask == expected_membership &&
                  local.stabilizer_indices == expected_indices &&
                  group->stabilizer_generators(local) == expected_generators,
                  "compact stabilizer membership changed the action matrices");
          require(local.canonical.automorphism_group_order == canonical.automorphism_group_order,
                  "explicit stabilizer disagrees with Bliss");
          for (const auto index : {std::size_t{0}, group->elements().size() - 1}) {
            std::vector<utsp::Mask> image;
            for (const auto row : rows) {
              image.push_back(utsp::linear_combination(row, group->elements()[index]));
            }
            require(group->classify(image).local_key == local.local_key,
                    "local orbit key changed under a support automorphism");
            require(utsp::canonicalize_marked_code(space, image).key_words == canonical.key_words,
                    "explicit orbit image changed the Bliss key");
            ++comparisons;
          }
        }
      }
      // An invariant subquotient makes complete differential traversal small.
      std::vector<utsp::Mask> invariant_basis;
      for (std::uint32_t bit = 0; bit < space.quotient_dimension(); ++bit) {
        auto proposed = invariant_basis;
        for (const auto& element : group->elements()) proposed.push_back(element[bit]);
        proposed = utsp::rref_basis(proposed);
        if (proposed.size() <= 5) invariant_basis = std::move(proposed);
        if (invariant_basis.size() == 5) break;
      }
      require(invariant_basis.size() >= 3, "fixture needs an invariant three-space");
      const auto restricted = restrict_space(space, invariant_basis);
      const auto restricted_root = utsp::canonicalize_marked_code(restricted, {});
      const auto restricted_group = utsp::SmallQuotientGroup::build(
          restricted_root, restricted.quotient_dimension());
      require(restricted_group.has_value(), "restricted small group did not close");
      incidence_checks += check_all_line_parent_orbits(
          *restricted_group, restricted.quotient_dimension());
      for (const auto target : {std::uint32_t{3}, std::uint32_t{5}}) {
        if (target > restricted.quotient_dimension()) continue;
        const std::array<std::uint64_t, 2> signatures = target == 3
            ? std::array<std::uint64_t, 2>{0x40, 0x41}
            : std::array<std::uint64_t, 2>{0x60000, 0x60001};
        const auto profile = utsp::build_tensor_restriction_profile(target, signatures, 1, 1);
        const auto late = utsp::enumerate_orbits_impl(
            restricted, target, 0, false, &profile, 3, 3, 1);
        const auto early = utsp::enumerate_orbits_impl(
            restricted, target, 0, false, &profile, 1, 3, 1,
            utsp::Q5ParentDomain::unrestricted, nullptr, true, false);
        const auto streamed = utsp::enumerate_orbits_impl(
            restricted, target, 0, false, &profile, 1, 3, 2,
            utsp::Q5ParentDomain::unrestricted, nullptr, true, true);
        require(late.canonical_keys == early.canonical_keys &&
                late.orbit_sizes == early.orbit_sizes &&
                late.weighted_subspace_count == early.weighted_subspace_count,
                "early target restrictions changed the exact orbit catalogue");
        require(early.canonical_keys == streamed.canonical_keys &&
                early.orbit_sizes == streamed.orbit_sizes &&
                early.weighted_subspace_count == streamed.weighted_subspace_count,
                "streaming two-row parents changed exact marked orbits");
        ++early_restriction_comparisons;
      }
      for (std::uint32_t q = 1; q <= 3; ++q) {
        for (const auto minimum : {std::uint32_t{0}, q}) {
          const auto reference = utsp::enumerate_orbits_impl(
              restricted, q, minimum, false, nullptr, 0, 0, 1,
              utsp::Q5ParentDomain::unrestricted, nullptr, false);
          const auto candidate = utsp::enumerate_orbits_impl(
              restricted, q, minimum, false, nullptr, 0, 0, 2,
              utsp::Q5ParentDomain::unrestricted, nullptr, true);
          require(reference.canonical_keys == candidate.canonical_keys &&
                  reference.orbit_sizes == candidate.orbit_sizes &&
                  reference.weighted_subspace_count == candidate.weighted_subspace_count,
                  "small-group traversal changed the exact orbit catalogue");
          ++enumerations;
        }
      }
    }
    std::cout << "{\"status\":\"pass\",\"image_comparisons\":" << comparisons
              << ",\"complete_enumeration_comparisons\":" << enumerations
              << ",\"canonical_line_incidence_checks\":" << incidence_checks
              << ",\"large_image_incidence_checks\":" << large_incidence_checks
              << ",\"image_order_1440_incidence_checks\":32"
              << ",\"early_restriction_comparisons\":" << early_restriction_comparisons << "}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
