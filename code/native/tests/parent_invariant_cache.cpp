// Include the implementation only in this test to check the private invariant.
#include "../src/marked_code_canonical.cpp"

#include <iostream>
#include <random>

int main() {
  try {
    std::vector<std::uint32_t> points;
    for (std::uint32_t point = 16; point < 64; ++point) points.push_back(point);
    const auto space = utsp::build_logical_label_space(points, 6, false);
    std::mt19937_64 random(20260906);
    std::uint64_t comparisons = 0;
    for (std::uint32_t dimension = 2; dimension <= 7; ++dimension) {
      for (std::uint32_t sample = 0; sample < 4; ++sample) {
        std::vector<utsp::Mask> children;
        while (children.size() < dimension) {
          const auto candidate = random() &
              ((utsp::Mask{1} << space.quotient_dimension()) - 1);
          auto proposed = children;
          proposed.push_back(candidate);
          if (utsp::gf2_rank(proposed) == proposed.size()) children = std::move(proposed);
        }
        std::vector<std::vector<utsp::Mask>> invariants;
        for (utsp::Mask functional = 1; functional < (utsp::Mask{1} << dimension);
             ++functional) {
          const std::array<utsp::Mask, 1> constraint{functional};
          const auto kernel = utsp::nullspace_basis(constraint, dimension);
          std::vector<utsp::Mask> parent;
          for (const auto word : kernel) {
            parent.push_back(utsp::linear_combination(word, children));
          }
          invariants.push_back(utsp::marked_code_invariant(space, parent));
        }
        auto ordered = invariants;
        std::sort(ordered.begin(), ordered.end());
        for (const auto index : {std::size_t{0}, ordered.size() / 2, ordered.size() - 1}) {
          const auto& parent = ordered[index];
          utsp::HyperplaneInvariantCache cache(space, children, parent);
          for (std::size_t i = 0; i < invariants.size(); ++i) {
            if (cache.smaller(i + 1) != (invariants[i] < parent)) {
              throw std::logic_error("cached invariant changed a parent comparison");
            }
            ++comparisons;
          }
          for (std::size_t i = invariants.size(); i > 0; --i) {
            if (cache.smaller(i) != (invariants[i - 1] < parent)) {
              throw std::logic_error("cached invariant depends on visitation order");
            }
            ++comparisons;
          }
        }
      }
    }
    std::cout << "{\"status\":\"pass\",\"hyperplane_comparisons\":"
              << comparisons << ",\"child_dimensions\":[2,7]}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
