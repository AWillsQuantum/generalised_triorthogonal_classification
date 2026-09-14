#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace utsp {

using Mask = std::uint64_t;

struct LabelSpace {
  std::vector<std::uint32_t> points;
  std::uint32_t ambient_dimension = 0;
  std::vector<Mask> coordinate_masks;
  std::vector<Mask> stabilizer_basis;
  std::vector<Mask> quotient_basis;
  std::vector<std::vector<Mask>> bilinear_form_rows;

  [[nodiscard]] std::uint32_t quotient_dimension() const {
    return static_cast<std::uint32_t>(quotient_basis.size());
  }
};

struct DistanceResult {
  std::optional<std::uint32_t> distance;
  std::uint64_t error_coefficient = 0;
};

struct IsotropicSummary {
  std::uint64_t subspace_count = 0;
  std::uint64_t tensor_checksum = 1469598103934665603ULL;
};

class DistanceContext {
 public:
  DistanceContext(
      std::span<const std::uint32_t> points,
      std::uint32_t ambient_dimension);

  [[nodiscard]] std::uint64_t weight_three_error_coefficient(
      std::span<const Mask> label_rows) const;

  [[nodiscard]] std::uint64_t weight_four_error_coefficient(
      std::span<const Mask> label_rows) const;

  [[nodiscard]] DistanceResult evaluate(
      std::span<const Mask> label_rows,
      std::uint32_t maximum_dp_dimension = 24) const;

 private:
  std::vector<std::uint32_t> points_;
  std::uint32_t ambient_dimension_ = 0;
  std::vector<Mask> weight_three_subsets_;
  mutable bool weight_four_ready_ = false;
  mutable std::vector<Mask> weight_four_subsets_;
};

[[nodiscard]] std::vector<Mask> rref_basis(std::span<const Mask> vectors);
[[nodiscard]] std::vector<Mask> nullspace_basis(
    std::span<const Mask> rows, std::uint32_t width);
[[nodiscard]] Mask reduce_vector(Mask vector, std::span<const Mask> basis);
[[nodiscard]] Mask linear_combination(
    Mask coefficients, std::span<const Mask> vectors);
[[nodiscard]] std::uint32_t gf2_rank(std::span<const Mask> vectors);

void validate_stabilizer_moments(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool require_even = false);

[[nodiscard]] LabelSpace build_logical_label_space(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool validate_moments = true,
    std::uint32_t minimum_distance = 3);

// Quotient-only construction for profiling and cache generation.  It omits
// coordinate-form reconstruction, which downstream cache consumers perform
// only for supports that survive their dimension gate.
[[nodiscard]] std::vector<Mask> build_logical_label_quotient_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool validate_moments = true,
    std::uint32_t minimum_distance = 3);

// Reconstruct the cheap, support-dependent parts of a label space from a
// previously certified quotient basis.  This avoids repeating the expensive
// low-weight-relation/nullspace calculation across logical dimensions.
[[nodiscard]] LabelSpace build_logical_label_space_from_quotient_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> quotient_basis,
    bool validate_moments = true);

[[nodiscard]] bool alternating_form_span_allows_isotropic_dimension(
    const LabelSpace& label_space, std::uint32_t dimension);

[[nodiscard]] std::vector<Mask> low_weight_relation_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t minimum_distance);

// Intersect a certified d=3 quotient with the additional low-weight error
// relations for a higher target distance.  Stabilizer rows already satisfy
// these relations, so the intersection can be computed in quotient
// coordinates without rebuilding the base label nullspace.
[[nodiscard]] LabelSpace restrict_logical_label_space_minimum_distance(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> distance_three_quotient_basis,
    std::uint32_t minimum_distance,
    bool validate_moments = true);

[[nodiscard]] std::vector<Mask>
restrict_logical_label_quotient_basis_minimum_distance(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> distance_three_quotient_basis,
    std::uint32_t minimum_distance,
    bool validate_moments = true);

using IsotropicVisitor =
    std::function<void(std::span<const Mask>, std::span<const Mask>)>;

std::uint64_t enumerate_totally_isotropic_subspaces(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    const IsotropicVisitor& visitor = {});

// Exact count-only path. Pivot sets are disjoint RREF cells, so they can be
// distributed independently without synchronization or duplicate subspaces.
std::uint64_t count_totally_isotropic_subspaces_parallel(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers);

[[nodiscard]] std::vector<std::uint64_t> cubic_tensor_words(
    std::span<const Mask> label_rows);

[[nodiscard]] std::uint64_t cubic_tensor_word(
    std::span<const Mask> label_rows);

[[nodiscard]] std::vector<Mask> tensor_radical_basis(
    std::span<const Mask> label_rows);

// The repeated-index part B(x,y)=T(x,x,y)=T(x,y,y) is the ordinary
// binary Gram form of the logical rows.
[[nodiscard]] std::vector<Mask> tensor_bilinear_form_rows(
    std::span<const Mask> label_rows);

[[nodiscard]] std::uint32_t tensor_bilinear_rank(
    std::span<const Mask> label_rows);

[[nodiscard]] bool tensor_bilinear_is_alternating(
    std::span<const Mask> label_rows);

// Count codimension-one logical subspaces on which the cubic tensor remains
// nondegenerate.  Intended for exact hereditary/primitive certification.
[[nodiscard]] std::uint64_t nondegenerate_tensor_hyperplane_count(
    std::span<const Mask> label_rows);

// Exact primitive test after the Gram-sector theorem: the tensor must be
// nondegenerate, have no nondegenerate hyperplane, and lie in one of the
// only Gram sectors in which a primitive tensor can occur.
[[nodiscard]] bool is_primitive_exceptional_tensor(
    std::span<const Mask> label_rows);

// Add the alternating form C=B+p tensor p, where p(x)=B(x,x).
// A subspace is isotropic for C exactly when B has rank at most one on it.
[[nodiscard]] LabelSpace rank_at_most_one_tensor_sector(
    const LabelSpace& label_space);

// Restrict the quotient to even-weight labels.  On this kernel B is
// alternating; full-rank B candidates form the second primitive sector.
[[nodiscard]] LabelSpace even_weight_tensor_sector(
    const LabelSpace& label_space);

[[nodiscard]] std::uint32_t intrinsic_tensor_dimension(
    std::span<const Mask> label_rows);

[[nodiscard]] DistanceResult z_distance_and_error_coefficient(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> label_rows,
    std::uint32_t maximum_dp_dimension = 24);

[[nodiscard]] IsotropicSummary summarize_isotropic_subspaces(
    const LabelSpace& label_space, std::uint32_t dimension);

[[nodiscard]] std::string hex_mask(Mask value);

}  // namespace utsp
