#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace utsp {

struct OriginTransporter {
  std::uint32_t origin = 0;
  std::vector<std::uint32_t> source_dual_basis;
  std::vector<std::uint32_t> image_dual_basis;
};

struct AffineOriginOrbit {
  std::uint32_t representative = 0;
  std::vector<OriginTransporter> members;
};

struct AffineOriginOrbitResult {
  std::uint32_t ambient_dimension = 0;
  std::vector<AffineOriginOrbit> orbits;
  std::uint64_t invariant_bucket_count = 0;
  std::uint64_t exact_equivalence_tests = 0;
  std::uint64_t exact_search_nodes = 0;
  std::uint64_t coarse_search_nodes = 0;
  std::uint64_t unmarked_refinement_search_nodes = 0;
  std::uint64_t refined_dual_search_nodes = 0;
  std::uint64_t primal_search_nodes = 0;
  std::uint64_t second_order_refinements = 0;
  bool used_refined_colors = false;
  bool used_marked_relations = false;
  bool used_primal_search = false;
  bool complete = true;
};

// Exact affine-stabilizer orbits on every point of F_2^m. The returned dual
// transporters map the orbit representative to each listed origin.
[[nodiscard]] AffineOriginOrbitResult affine_origin_orbits(
    std::span<const std::uint32_t> support_points,
    std::uint32_t ambient_dimension,
    std::uint64_t node_limit = 0);

}  // namespace utsp
