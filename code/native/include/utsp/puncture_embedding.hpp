#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace utsp {

struct HereditaryPunctureSearchResult {
  std::vector<std::uint64_t> puncture_index_masks;
  std::uint64_t search_nodes = 0;
  std::uint64_t profile_rejections = 0;
  std::uint64_t complete_embeddings = 0;
  std::uint64_t unique_direction_spaces = 0;
  std::uint64_t unique_embedded_patterns = 0;
  std::uint64_t duplicate_embedded_patterns = 0;
  std::uint64_t candidate_cosets = 0;
  bool histogram_rejected = false;
  bool complete = true;
};

[[nodiscard]] std::vector<std::uint32_t> pair_difference_profile(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension);

// Find every support subset that is an affine image of template_points and
// whose direction-space cosets outside the distinguished flat are singletons.
[[nodiscard]] HereditaryPunctureSearchResult find_hereditary_punctures(
    std::span<const std::uint32_t> support_points,
    std::uint32_t ambient_dimension,
    std::span<const std::uint32_t> template_points,
    std::uint32_t template_dimension,
    std::uint64_t node_limit = 0);

}  // namespace utsp
