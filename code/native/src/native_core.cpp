#include "utsp/native_core.hpp"

#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <atomic>
#include <array>
#include <bit>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_set>
#include <unordered_map>
#include <utility>

namespace utsp {
namespace {

constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

[[nodiscard]] std::uint32_t pivot_of(Mask value) {
  if (value == 0) throw std::invalid_argument("zero vector has no pivot");
  return 63U - static_cast<std::uint32_t>(std::countl_zero(value));
}

[[nodiscard]] bool dot(Mask left, Mask right) {
  return (std::popcount(left & right) & 1U) != 0;
}

[[nodiscard]] Mask width_mask(std::uint32_t width) {
  if (width > 64) throw std::invalid_argument("GF(2) width exceeds 64");
  return width == 64 ? std::numeric_limits<Mask>::max()
                     : ((Mask{1} << width) - 1);
}

[[nodiscard]] std::uint32_t row_rank(std::span<const Mask> rows) {
  std::array<Mask, 64> basis{};
  std::uint32_t rank = 0;
  for (Mask value : rows) {
    while (value != 0) {
      const auto pivot = pivot_of(value);
      if (basis[pivot] != 0) {
        value ^= basis[pivot];
        continue;
      }
      basis[pivot] = value;
      ++rank;
      break;
    }
  }
  return rank;
}

class OrthogonalityCache {
 public:
  explicit OrthogonalityCache(const LabelSpace& label_space)
      : forms_(label_space.bilinear_form_rows),
        width_(label_space.quotient_dimension()) {
    if (width_ <= 21) {
      constraints_.resize(std::size_t{1} << width_);
      ready_.assign(std::size_t{1} << width_, false);
    }
  }

  [[nodiscard]] bool orthogonal(Mask left, Mask right) {
    if (forms_.empty()) return true;
    const auto& constraints = constraints_for(left);
    return std::none_of(
        constraints.begin(), constraints.end(),
        [right](Mask constraint) { return dot(constraint, right); });
  }

  [[nodiscard]] std::uint32_t common_orthogonal_dimension(
      std::span<const Mask> rows) {
    return width_ - static_cast<std::uint32_t>(
                        common_constraints(rows).size());
  }

  [[nodiscard]] std::vector<Mask> common_constraints(
      std::span<const Mask> rows) {
    std::vector<Mask> constraints;
    constraints.reserve(rows.size() * forms_.size());
    for (const auto row : rows) {
      const auto& row_constraints = constraints_for(row);
      constraints.insert(
          constraints.end(), row_constraints.begin(), row_constraints.end());
    }
    return rref_basis(constraints);
  }

  [[nodiscard]] const std::vector<Mask>& constraints_for_vector(Mask row) {
    return constraints_for(row);
  }

  [[nodiscard]] bool residual_form_span_allows_target(
      std::span<const Mask> common_constraints,
      std::uint32_t selected_dimension,
      std::uint32_t target_dimension) const {
    const auto common_dimension =
        width_ - static_cast<std::uint32_t>(common_constraints.size());
    if (common_dimension < target_dimension) return false;
    const auto quotient_dimension = common_dimension - selected_dimension;
    const auto remaining_dimension = target_dimension - selected_dimension;
    const auto required_rank =
        2U * (quotient_dimension - remaining_dimension + 1U);
    if (required_rank > 2U * (quotient_dimension / 2U) || forms_.empty()) {
      return true;
    }
    if (forms_.size() >= 63 ||
        width_ + common_constraints.size() > 64) {
      throw std::overflow_error("residual alternating pencil exceeds uint64_t");
    }

    const auto restricted_rank = [&](std::span<const Mask> form) {
      std::vector<Mask> bordered;
      bordered.reserve(width_ + common_constraints.size());
      for (std::uint32_t row = 0; row < width_; ++row) {
        Mask value = form[row];
        for (std::size_t constraint = 0;
             constraint < common_constraints.size(); ++constraint) {
          if ((common_constraints[constraint] >> row) & Mask{1}) {
            value |= Mask{1} << (width_ + constraint);
          }
        }
        bordered.push_back(value);
      }
      bordered.insert(
          bordered.end(), common_constraints.begin(), common_constraints.end());
      const auto bordered_rank = row_rank(bordered);
      const auto constraint_rank = static_cast<std::uint32_t>(
          common_constraints.size());
      if (bordered_rank < 2U * constraint_rank) {
        throw std::logic_error("bordered restriction rank is inconsistent");
      }
      return bordered_rank - 2U * constraint_rank;
    };

    std::vector<Mask> combination(width_, 0);
    std::uint64_t previous_gray = 0;
    const auto combination_count = std::uint64_t{1} << forms_.size();
    for (std::uint64_t index = 1; index < combination_count; ++index) {
      const auto gray = index ^ (index >> 1U);
      const auto changed = gray ^ previous_gray;
      const auto form_index = static_cast<std::size_t>(
          std::countr_zero(changed));
      for (std::uint32_t row = 0; row < width_; ++row) {
        combination[row] ^= forms_[form_index][row];
      }
      previous_gray = gray;
      if (restricted_rank(combination) >= required_rank) return false;
    }
    return true;
  }

 private:
  [[nodiscard]] const std::vector<Mask>& constraints_for(Mask left) {
    if (!constraints_.empty()) {
      const auto index = static_cast<std::size_t>(left);
      if (!ready_[index]) {
        constraints_[index] = build_constraints(left);
        ready_[index] = true;
      }
      return constraints_[index];
    }
    fallback_ = build_constraints(left);
    return fallback_;
  }

  [[nodiscard]] std::vector<Mask> build_constraints(Mask left) const {
    std::vector<Mask> contracted;
    contracted.reserve(forms_.size());
    for (const auto& form : forms_) {
      contracted.push_back(linear_combination(left, form));
    }
    return rref_basis(contracted);
  }

  const std::vector<std::vector<Mask>>& forms_;
  std::uint32_t width_;
  std::vector<std::vector<Mask>> constraints_;
  std::vector<bool> ready_;
  std::vector<Mask> fallback_;
};

void enumerate_pivot_sets(
    std::uint32_t width,
    std::uint32_t dimension,
    std::uint32_t next,
    std::vector<std::uint32_t>& pivots,
    const std::function<void(std::span<const std::uint32_t>)>& visitor) {
  if (pivots.size() == dimension) {
    visitor(pivots);
    return;
  }
  const auto remaining = dimension - static_cast<std::uint32_t>(pivots.size());
  for (std::uint32_t pivot = next; pivot + remaining <= width; ++pivot) {
    pivots.push_back(pivot);
    enumerate_pivot_sets(width, dimension, pivot + 1, pivots, visitor);
    pivots.pop_back();
  }
}

struct AffineRowSpace {
  bool consistent = true;
  Mask particular = 0;
  std::vector<Mask> nullspace_basis;
};

[[nodiscard]] AffineRowSpace constrained_row_space(
    std::uint32_t pivot,
    std::span<const std::uint32_t> pivots,
    std::uint32_t width,
    std::span<const Mask> constraints) {
  Mask pivot_mask = 0;
  for (const auto value : pivots) pivot_mask |= Mask{1} << value;
  Mask free_mask = 0;
  for (std::uint32_t column = pivot + 1; column < width; ++column) {
    if (((pivot_mask >> column) & 1U) == 0) {
      free_mask |= Mask{1} << column;
    }
  }

  std::array<Mask, 64> equations{};
  std::array<bool, 64> right_hand_sides{};
  for (const auto constraint : constraints) {
    Mask equation = constraint & free_mask;
    bool right_hand_side =
        ((constraint >> pivot) & Mask{1}) != 0;
    while (equation != 0) {
      const auto equation_pivot = pivot_of(equation);
      if (equations[equation_pivot] != 0) {
        equation ^= equations[equation_pivot];
        right_hand_side =
            right_hand_side != right_hand_sides[equation_pivot];
        continue;
      }
      equations[equation_pivot] = equation;
      right_hand_sides[equation_pivot] = right_hand_side;
      for (std::uint32_t other = 0; other < width; ++other) {
        if (other == equation_pivot || equations[other] == 0 ||
            ((equations[other] >> equation_pivot) & Mask{1}) == 0) {
          continue;
        }
        equations[other] ^= equation;
        right_hand_sides[other] =
            right_hand_sides[other] != right_hand_side;
      }
      break;
    }
    if (equation == 0 && right_hand_side) {
      return AffineRowSpace{false, 0, {}};
    }
  }

  AffineRowSpace result;
  result.particular = Mask{1} << pivot;
  Mask equation_pivots = 0;
  for (std::uint32_t column = 0; column < width; ++column) {
    if (equations[column] == 0) continue;
    equation_pivots |= Mask{1} << column;
    if (right_hand_sides[column]) result.particular |= Mask{1} << column;
  }
  Mask free_variables = free_mask & ~equation_pivots;
  while (free_variables != 0) {
    const auto free = static_cast<std::uint32_t>(
        std::countr_zero(free_variables));
    Mask basis_vector = Mask{1} << free;
    for (std::uint32_t column = 0; column < width; ++column) {
      if (equations[column] != 0 &&
          ((equations[column] >> free) & Mask{1}) != 0) {
        basis_vector |= Mask{1} << column;
      }
    }
    result.nullspace_basis.push_back(basis_vector);
    free_variables &= free_variables - 1;
  }
  return result;
}

void enumerate_rows(
    std::size_t index,
    std::span<const std::uint32_t> pivots,
    std::uint32_t width,
    std::vector<Mask>& rows,
    OrthogonalityCache& cache,
    std::span<const Mask> common_constraints,
    bool residual_rank_pruning,
    const std::function<void(std::span<const Mask>)>& visitor) {
  if (index != 0 &&
      width - static_cast<std::uint32_t>(common_constraints.size()) <
          pivots.size()) {
    return;
  }
  if (index == pivots.size()) {
    visitor(rows);
    return;
  }
  if (residual_rank_pruning && index >= 2 &&
      !cache.residual_form_span_allows_target(
          common_constraints,
          static_cast<std::uint32_t>(index),
          static_cast<std::uint32_t>(pivots.size()))) {
    return;
  }
  const auto row_space = constrained_row_space(
      pivots[index], pivots, width, common_constraints);
  if (!row_space.consistent) return;
  if (row_space.nullspace_basis.size() >= 63) {
    throw std::overflow_error("RREF row-choice count exceeds uint64_t");
  }
  const auto choice_count =
      std::uint64_t{1} << row_space.nullspace_basis.size();
  std::vector<Mask> admissible_rows;
  admissible_rows.reserve(static_cast<std::size_t>(choice_count));
  for (std::uint64_t choice = 0; choice < choice_count; ++choice) {
    admissible_rows.push_back(
        row_space.particular ^
        linear_combination(choice, row_space.nullspace_basis));
  }
  std::sort(admissible_rows.begin(), admissible_rows.end());
  for (const Mask row : admissible_rows) {
    rows[index] = row;
    std::vector<Mask> child_constraints(
        common_constraints.begin(), common_constraints.end());
    const auto& row_constraints = cache.constraints_for_vector(row);
    child_constraints.insert(
        child_constraints.end(), row_constraints.begin(), row_constraints.end());
    child_constraints = rref_basis(child_constraints);
    enumerate_rows(
        index + 1,
        pivots,
        width,
        rows,
        cache,
        child_constraints,
        residual_rank_pruning,
        visitor);
  }
}

void set_tensor_bit(std::vector<std::uint64_t>& words, std::size_t bit) {
  words[bit / 64] |= std::uint64_t{1} << (bit % 64);
}

void fnv_word(std::uint64_t& state, std::uint64_t word) {
  for (unsigned byte = 0; byte < 8; ++byte) {
    state ^= (word >> (8 * byte)) & 0xffU;
    state *= kFnvPrime;
  }
}

[[nodiscard]] LabelSpace initialize_label_space(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool validate_moments) {
  if (points.empty()) throw std::invalid_argument("support is empty");
  if (points.size() > 64) {
    throw std::invalid_argument("support length exceeds 64-bit row masks");
  }
  if (ambient_dimension >= 32) {
    throw std::invalid_argument("point representation supports h < 32");
  }
  std::unordered_set<std::uint32_t> distinct;
  std::vector<Mask> point_vectors;
  point_vectors.reserve(points.size());
  for (const auto point : points) {
    if (point >= (std::uint64_t{1} << ambient_dimension)) {
      throw std::invalid_argument("point is outside ambient space");
    }
    if (!distinct.insert(point).second) {
      throw std::invalid_argument("support contains a repeated point");
    }
    point_vectors.push_back(point);
  }
  if (gf2_rank(point_vectors) != ambient_dimension) {
    throw std::invalid_argument("support does not span ambient space");
  }
  if (validate_moments) {
    validate_stabilizer_moments(points, ambient_dimension);
  }

  LabelSpace result;
  result.points.assign(points.begin(), points.end());
  result.ambient_dimension = ambient_dimension;
  result.coordinate_masks.reserve(ambient_dimension);
  for (std::uint32_t coordinate = 0; coordinate < ambient_dimension;
       ++coordinate) {
    Mask mask = 0;
    for (std::size_t index = 0; index < points.size(); ++index) {
      if ((points[index] >> coordinate) & 1U) mask |= Mask{1} << index;
    }
    result.coordinate_masks.push_back(mask);
  }
  result.stabilizer_basis = rref_basis(result.coordinate_masks);
  if (result.stabilizer_basis.size() != ambient_dimension) {
    throw std::logic_error("coordinate masks unexpectedly lost rank");
  }
  return result;
}

[[nodiscard]] std::vector<Mask> base_label_constraints(
    const LabelSpace& label_space) {
  std::vector<Mask> constraints = label_space.coordinate_masks;
  for (std::uint32_t left = 0; left < label_space.ambient_dimension; ++left) {
    for (std::uint32_t right = left + 1;
         right < label_space.ambient_dimension;
         ++right) {
      constraints.push_back(
          label_space.coordinate_masks[left] &
          label_space.coordinate_masks[right]);
    }
  }
  const auto zero =
      std::find(label_space.points.begin(), label_space.points.end(), 0U);
  if (zero != label_space.points.end()) {
    constraints.push_back(
        Mask{1} << std::distance(label_space.points.begin(), zero));
  }
  return constraints;
}

void populate_bilinear_forms(LabelSpace& label_space) {
  label_space.bilinear_form_rows.clear();
  label_space.bilinear_form_rows.reserve(label_space.ambient_dimension);
  for (const Mask coordinate_mask : label_space.coordinate_masks) {
    std::vector<Mask> form;
    form.reserve(label_space.quotient_basis.size());
    for (const Mask left : label_space.quotient_basis) {
      Mask row = 0;
      for (std::size_t right_index = 0;
           right_index < label_space.quotient_basis.size();
           ++right_index) {
        if (dot(
                left & label_space.quotient_basis[right_index],
                coordinate_mask)) {
          row |= Mask{1} << right_index;
        }
      }
      form.push_back(row);
    }
    for (std::size_t index = 0; index < form.size(); ++index) {
      if ((form[index] >> index) & 1U) {
        throw std::logic_error("induced bilinear form is not alternating");
      }
    }
    label_space.bilinear_form_rows.push_back(std::move(form));
  }
}

void populate_label_quotient_basis(
    LabelSpace& label_space,
    std::uint32_t minimum_distance) {
  auto constraints = base_label_constraints(label_space);
  if (minimum_distance > 3) {
    const auto relations =
        low_weight_relation_basis(label_space.points, minimum_distance);
    constraints.insert(constraints.end(), relations.begin(), relations.end());
  }

  const auto nullspace = nullspace_basis(constraints, label_space.points.size());
  for (const Mask constraint : constraints) {
    for (const Mask row : label_space.stabilizer_basis) {
      if (dot(constraint, row)) {
        throw std::logic_error("stabilizer space is not contained in label space");
      }
    }
  }

  std::vector<Mask> span_basis = label_space.stabilizer_basis;
  for (const Mask row : nullspace) {
    const Mask residual = reduce_vector(row, span_basis);
    if (residual == 0) continue;
    label_space.quotient_basis.push_back(residual);
    span_basis.push_back(residual);
    span_basis = rref_basis(span_basis);
  }
  if (label_space.quotient_basis.size() + label_space.ambient_dimension !=
      nullspace.size()) {
    throw std::logic_error("failed to construct label-space quotient");
  }
}

void install_label_quotient_basis(
    LabelSpace& label_space,
    std::span<const Mask> quotient_basis) {
  const auto constraints = base_label_constraints(label_space);
  const auto support_mask =
      width_mask(static_cast<std::uint32_t>(label_space.points.size()));
  std::vector<Mask> span_basis = label_space.stabilizer_basis;
  label_space.quotient_basis.reserve(quotient_basis.size());
  for (const Mask row : quotient_basis) {
    if ((row & ~support_mask) != 0) {
      throw std::invalid_argument("cached quotient row exceeds support width");
    }
    for (const Mask constraint : constraints) {
      if (dot(constraint, row)) {
        throw std::invalid_argument(
            "cached quotient row violates a base label constraint");
      }
    }
    if (reduce_vector(row, span_basis) == 0) {
      throw std::invalid_argument(
          "cached quotient basis is dependent modulo the stabilizer");
    }
    label_space.quotient_basis.push_back(row);
    span_basis.push_back(row);
    span_basis = rref_basis(span_basis);
  }
  if (label_space.quotient_basis.size() + label_space.ambient_dimension >
      label_space.points.size()) {
    throw std::invalid_argument("cached quotient dimension is impossible");
  }
}

}  // namespace

std::vector<Mask> rref_basis(std::span<const Mask> vectors) {
  std::map<std::uint32_t, Mask, std::greater<>> basis;
  for (const Mask original : vectors) {
    Mask value = original;
    for (const auto& [pivot, row] : basis) {
      if ((value >> pivot) & 1U) value ^= row;
    }
    if (value == 0) continue;
    const auto pivot = pivot_of(value);
    for (auto& [other_pivot, row] : basis) {
      (void)other_pivot;
      if ((row >> pivot) & 1U) row ^= value;
    }
    basis.emplace(pivot, value);
  }

  std::vector<Mask> result;
  result.reserve(basis.size());
  for (const auto& [pivot, row] : basis) {
    (void)pivot;
    result.push_back(row);
  }
  return result;
}

std::vector<Mask> nullspace_basis(
    std::span<const Mask> rows, std::uint32_t width) {
  const Mask mask = width_mask(width);
  std::vector<Mask> normalized;
  normalized.reserve(rows.size());
  for (const Mask row : rows) normalized.push_back(row & mask);
  const auto reduced = rref_basis(normalized);

  Mask pivot_mask = 0;
  for (const Mask row : reduced) pivot_mask |= Mask{1} << pivot_of(row);

  std::vector<Mask> basis;
  for (std::uint32_t free = 0; free < width; ++free) {
    if ((pivot_mask >> free) & 1U) continue;
    Mask vector = Mask{1} << free;
    for (const Mask row : reduced) {
      if (dot(row, vector)) vector |= Mask{1} << pivot_of(row);
    }
    basis.push_back(vector);
  }
  return basis;
}

Mask reduce_vector(Mask vector, std::span<const Mask> basis) {
  auto reduced_basis = rref_basis(basis);
  for (const Mask row : reduced_basis) {
    const auto pivot = pivot_of(row);
    if ((vector >> pivot) & 1U) vector ^= row;
  }
  return vector;
}

Mask linear_combination(Mask coefficients, std::span<const Mask> vectors) {
  Mask result = 0;
  while (coefficients != 0) {
    const auto index = static_cast<std::size_t>(std::countr_zero(coefficients));
    if (index >= vectors.size()) {
      throw std::invalid_argument("coefficient lies outside vector list");
    }
    result ^= vectors[index];
    coefficients &= coefficients - 1;
  }
  return result;
}

std::uint32_t gf2_rank(std::span<const Mask> vectors) {
  return static_cast<std::uint32_t>(rref_basis(vectors).size());
}

void validate_stabilizer_moments(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool require_even) {
  if (ambient_dimension >= 32) {
    throw std::invalid_argument("point representation supports h < 32");
  }
  if (require_even && (points.size() & 1U)) {
    throw std::invalid_argument("unital support has odd cardinality");
  }
  for (std::uint32_t first = 0; first < ambient_dimension; ++first) {
    unsigned parity = 0;
    for (const auto point : points) parity ^= (point >> first) & 1U;
    if (parity) throw std::invalid_argument("degree-one moment is odd");
    for (std::uint32_t second = first + 1; second < ambient_dimension; ++second) {
      parity = 0;
      for (const auto point : points) {
        parity ^= ((point >> first) & 1U) & ((point >> second) & 1U);
      }
      if (parity) throw std::invalid_argument("degree-two moment is odd");
      for (std::uint32_t third = second + 1; third < ambient_dimension; ++third) {
        parity = 0;
        for (const auto point : points) {
          parity ^= ((point >> first) & 1U) & ((point >> second) & 1U) &
                    ((point >> third) & 1U);
        }
        if (parity) throw std::invalid_argument("degree-three moment is odd");
      }
    }
  }
}

LabelSpace build_logical_label_space(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool validate_moments,
    std::uint32_t minimum_distance) {
  if (minimum_distance < 3 || minimum_distance > 5) {
    throw std::invalid_argument("implemented distance filtration is d=3,4,5");
  }

  auto result = initialize_label_space(
      points, ambient_dimension, validate_moments);
  populate_label_quotient_basis(result, minimum_distance);
  populate_bilinear_forms(result);
  return result;
}

std::vector<Mask> build_logical_label_quotient_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    bool validate_moments,
    std::uint32_t minimum_distance) {
  if (minimum_distance < 3 || minimum_distance > 5) {
    throw std::invalid_argument("implemented distance filtration is d=3,4,5");
  }
  auto label_space = initialize_label_space(
      points, ambient_dimension, validate_moments);
  populate_label_quotient_basis(label_space, minimum_distance);
  return std::move(label_space.quotient_basis);
}

LabelSpace build_logical_label_space_from_quotient_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> quotient_basis,
    bool validate_moments) {
  auto result = initialize_label_space(
      points, ambient_dimension, validate_moments);
  install_label_quotient_basis(result, quotient_basis);
  populate_bilinear_forms(result);
  return result;
}

bool alternating_form_span_allows_isotropic_dimension(
    const LabelSpace& label_space, std::uint32_t dimension) {
  const auto width = label_space.quotient_dimension();
  if (dimension == 0 || dimension > width) {
    throw std::invalid_argument("logical dimension is outside label quotient");
  }

  // An alternating form of rank r admits isotropic subspaces of dimension at
  // most width-r/2.  Every common isotropic subspace is isotropic for every
  // form in the linear span, so one sufficiently high-rank combination rules
  // out the requested dimension.
  const auto required_rank = 2U * (width - dimension + 1U);
  const auto maximum_alternating_rank = 2U * (width / 2U);
  if (required_rank > maximum_alternating_rank ||
      label_space.bilinear_form_rows.empty()) {
    return true;
  }
  if (label_space.bilinear_form_rows.size() >= 63) {
    throw std::overflow_error("bilinear-form pencil exceeds uint64_t");
  }
  for (const auto& form : label_space.bilinear_form_rows) {
    if (form.size() != width) {
      throw std::invalid_argument("bilinear form has the wrong dimension");
    }
    if (row_rank(form) >= required_rank) return false;
  }

  std::vector<Mask> combination(width, 0);
  const auto combination_count =
      std::uint64_t{1} << label_space.bilinear_form_rows.size();
  std::uint64_t previous_gray = 0;
  for (std::uint64_t index = 1; index < combination_count; ++index) {
    const auto gray = index ^ (index >> 1U);
    const auto changed = gray ^ previous_gray;
    const auto form_index = static_cast<std::size_t>(std::countr_zero(changed));
    const auto& form = label_space.bilinear_form_rows[form_index];
    for (std::size_t row = 0; row < combination.size(); ++row) {
      combination[row] ^= form[row];
    }
    previous_gray = gray;
    if (std::has_single_bit(gray)) continue;
    if (row_rank(combination) >= required_rank) return false;
  }
  return true;
}

std::vector<Mask> low_weight_relation_basis(
    std::span<const std::uint32_t> points,
    std::uint32_t minimum_distance) {
  if (points.size() > 64) {
    throw std::invalid_argument("relation basis supports at most 64 points");
  }
  if (minimum_distance < 1 || minimum_distance > 5) {
    throw std::invalid_argument("implemented relation weights stop below five");
  }
  std::vector<Mask> relations;
  if (minimum_distance > 1) {
    for (std::size_t index = 0; index < points.size(); ++index) {
      if (points[index] == 0) relations.push_back(Mask{1} << index);
    }
  }
  if (minimum_distance > 2) {
    std::unordered_map<std::uint32_t, std::size_t> first_index;
    for (std::size_t index = 0; index < points.size(); ++index) {
      const auto [found, inserted] = first_index.emplace(points[index], index);
      if (!inserted) {
        relations.push_back((Mask{1} << found->second) | (Mask{1} << index));
      }
    }
  }
  if (minimum_distance > 3) {
    std::unordered_map<std::uint32_t, std::size_t> point_index;
    for (std::size_t index = 0; index < points.size(); ++index) {
      point_index.emplace(points[index], index);
    }
    for (std::size_t left = 0; left < points.size(); ++left) {
      for (std::size_t right = left + 1; right < points.size(); ++right) {
        const auto found = point_index.find(points[left] ^ points[right]);
        if (found == point_index.end() || found->second <= right) continue;
        relations.push_back(
            (Mask{1} << left) | (Mask{1} << right) |
            (Mask{1} << found->second));
      }
    }
  }
  if (minimum_distance > 4) {
    std::unordered_map<std::uint32_t, std::vector<Mask>> pairs_by_xor;
    for (std::size_t left = 0; left < points.size(); ++left) {
      for (std::size_t right = left + 1; right < points.size(); ++right) {
        pairs_by_xor[points[left] ^ points[right]].push_back(
            (Mask{1} << left) | (Mask{1} << right));
      }
    }
    for (const auto& [syndrome, pairs] : pairs_by_xor) {
      (void)syndrome;
      for (std::size_t left = 0; left < pairs.size(); ++left) {
        for (std::size_t right = left + 1; right < pairs.size(); ++right) {
          if ((pairs[left] & pairs[right]) == 0) {
            relations.push_back(pairs[left] | pairs[right]);
          }
        }
      }
    }
  }
  return rref_basis(relations);
}

LabelSpace restrict_logical_label_space_minimum_distance(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> distance_three_quotient_basis,
    std::uint32_t minimum_distance,
    bool validate_moments) {
  const auto restricted_basis =
      restrict_logical_label_quotient_basis_minimum_distance(
          points,
          ambient_dimension,
          distance_three_quotient_basis,
          minimum_distance,
          validate_moments);
  return build_logical_label_space_from_quotient_basis(
      points,
      ambient_dimension,
      restricted_basis,
      false);
}

std::vector<Mask> restrict_logical_label_quotient_basis_minimum_distance(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> distance_three_quotient_basis,
    std::uint32_t minimum_distance,
    bool validate_moments) {
  if (minimum_distance < 4 || minimum_distance > 5) {
    throw std::invalid_argument(
        "derived label-space filtration supports target distance 4 or 5");
  }
  if (distance_three_quotient_basis.size() > 64) {
    throw std::invalid_argument("parent label quotient exceeds 64 dimensions");
  }

  auto parent = initialize_label_space(
      points, ambient_dimension, validate_moments);
  install_label_quotient_basis(parent, distance_three_quotient_basis);
  const auto relations = low_weight_relation_basis(points, minimum_distance);
  std::vector<Mask> quotient_constraints;
  quotient_constraints.reserve(relations.size());
  for (const Mask relation : relations) {
    for (const Mask stabilizer : parent.stabilizer_basis) {
      if (dot(relation, stabilizer)) {
        throw std::logic_error(
            "low-weight relation does not annihilate the stabilizer");
      }
    }
    Mask constraint = 0;
    for (std::size_t index = 0;
         index < parent.quotient_basis.size();
         ++index) {
      if (dot(relation, parent.quotient_basis[index])) {
        constraint |= Mask{1} << index;
      }
    }
    if (constraint != 0) quotient_constraints.push_back(constraint);
  }

  const auto quotient_kernel = nullspace_basis(
      quotient_constraints,
      static_cast<std::uint32_t>(parent.quotient_basis.size()));
  std::vector<Mask> restricted_basis;
  restricted_basis.reserve(quotient_kernel.size());
  for (const Mask coefficients : quotient_kernel) {
    restricted_basis.push_back(
        linear_combination(coefficients, parent.quotient_basis));
  }
  auto result = initialize_label_space(points, ambient_dimension, false);
  install_label_quotient_basis(result, restricted_basis);
  for (const Mask relation : relations) {
    for (const Mask row : result.quotient_basis) {
      if (dot(relation, row)) {
        throw std::logic_error(
            "derived label quotient violates a low-weight relation");
      }
    }
  }
  return std::move(result.quotient_basis);
}

std::uint64_t enumerate_totally_isotropic_subspaces(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    const IsotropicVisitor& visitor) {
  const auto width = label_space.quotient_dimension();
  if (dimension == 0 || dimension > width) {
    throw std::invalid_argument("logical dimension is outside label quotient");
  }
  if (!alternating_form_span_allows_isotropic_dimension(
          label_space, dimension)) {
    return 0;
  }

  OrthogonalityCache cache(label_space);
  std::uint64_t count = 0;
  std::vector<std::uint32_t> pivots;
  pivots.reserve(dimension);
  enumerate_pivot_sets(
      width, dimension, 0, pivots,
      [&](std::span<const std::uint32_t> selected_pivots) {
        std::vector<Mask> rows(dimension);
        enumerate_rows(
            0, selected_pivots, width, rows, cache, {}, false,
            [&](std::span<const Mask> ascending_rows) {
              ++count;
              if (!visitor) return;
              const auto quotient_rows = rref_basis(ascending_rows);
              std::vector<Mask> label_rows;
              label_rows.reserve(dimension);
              for (const Mask row : quotient_rows) {
                label_rows.push_back(
                    linear_combination(row, label_space.quotient_basis));
              }
              visitor(quotient_rows, label_rows);
            });
      });
  return count;
}

std::uint64_t count_totally_isotropic_subspaces_parallel(
    const LabelSpace& label_space,
    std::uint32_t dimension,
    std::uint32_t workers) {
  const auto width = label_space.quotient_dimension();
  if (dimension == 0 || dimension > width) {
    throw std::invalid_argument("logical dimension is outside label quotient");
  }
  if (workers == 0) {
    throw std::invalid_argument("isotropic counter worker count must be positive");
  }
  if (workers == 1) {
    return enumerate_totally_isotropic_subspaces(label_space, dimension);
  }
  if (!alternating_form_span_allows_isotropic_dimension(
          label_space, dimension)) {
    return 0;
  }

  std::vector<std::vector<std::uint32_t>> pivot_sets;
  std::vector<std::uint32_t> pivots;
  pivots.reserve(dimension);
  enumerate_pivot_sets(
      width,
      dimension,
      0,
      pivots,
      [&](std::span<const std::uint32_t> selected) {
        pivot_sets.emplace_back(selected.begin(), selected.end());
      });
  const auto worker_count = static_cast<std::size_t>(
      std::min<std::uint64_t>(workers, pivot_sets.size()));
  std::vector<std::uint64_t> counts(worker_count, 0);
  std::vector<std::exception_ptr> errors(worker_count);
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> threads;
  threads.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    threads.emplace_back([&, worker] {
      try {
        OrthogonalityCache cache(label_space);
        std::vector<Mask> rows(dimension);
        while (true) {
          const auto index = next.fetch_add(1);
          if (index >= pivot_sets.size()) return;
          enumerate_rows(
              0,
              pivot_sets[index],
              width,
              rows,
              cache,
              {},
              true,
              [&](std::span<const Mask>) { ++counts[worker]; });
        }
      } catch (...) {
        errors[worker] = std::current_exception();
        next.store(pivot_sets.size());
      }
    });
  }
  for (auto& thread : threads) thread.join();
  for (const auto& error : errors) {
    if (error) std::rethrow_exception(error);
  }
  return std::accumulate(counts.begin(), counts.end(), std::uint64_t{0});
}

std::vector<std::uint64_t> cubic_tensor_words(
    std::span<const Mask> label_rows) {
  const std::size_t q = label_rows.size();
  const std::size_t terms = q + q * (q - 1) / 2 + q * (q - 1) * (q - 2) / 6;
  std::vector<std::uint64_t> words((terms + 63) / 64, 0);
  std::size_t bit = 0;
  for (std::size_t first = 0; first < q; ++first, ++bit) {
    if (std::popcount(label_rows[first]) & 1U) set_tensor_bit(words, bit);
  }
  for (std::size_t first = 0; first < q; ++first) {
    for (std::size_t second = first + 1; second < q; ++second, ++bit) {
      if (std::popcount(label_rows[first] & label_rows[second]) & 1U) {
        set_tensor_bit(words, bit);
      }
    }
  }
  for (std::size_t first = 0; first < q; ++first) {
    for (std::size_t second = first + 1; second < q; ++second) {
      for (std::size_t third = second + 1; third < q; ++third, ++bit) {
        if (std::popcount(
                label_rows[first] & label_rows[second] & label_rows[third]) &
            1U) {
          set_tensor_bit(words, bit);
        }
      }
    }
  }
  return words;
}

std::uint64_t cubic_tensor_word(std::span<const Mask> label_rows) {
  const std::size_t q = label_rows.size();
  const std::size_t terms = q + q * (q - 1) / 2 + q * (q - 1) * (q - 2) / 6;
  if (terms > 64) {
    throw std::invalid_argument("cubic tensor does not fit in one word");
  }
  const auto words = cubic_tensor_words(label_rows);
  return words.empty() ? 0 : words[0];
}

std::vector<Mask> tensor_radical_basis(std::span<const Mask> label_rows) {
  const auto q = static_cast<std::uint32_t>(label_rows.size());
  if (q > 64) throw std::invalid_argument("tensor radical supports q <= 64");
  if (q >= 1 && q <= 7) {
    const auto packed = cubic_tensor_radical_basis_from_signature(
        q, cubic_tensor_word(label_rows));
    return {packed.begin(), packed.end()};
  }
  std::vector<Mask> constraints;
  constraints.reserve(static_cast<std::size_t>(q) * (q + 1) / 2);
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first; second < q; ++second) {
      Mask constraint = 0;
      const Mask overlap = label_rows[first] & label_rows[second];
      for (std::uint32_t third = 0; third < q; ++third) {
        if (std::popcount(overlap & label_rows[third]) & 1U) {
          constraint |= Mask{1} << third;
        }
      }
      constraints.push_back(constraint);
    }
  }
  return nullspace_basis(constraints, q);
}

std::vector<Mask> tensor_bilinear_form_rows(
    std::span<const Mask> label_rows) {
  if (label_rows.size() > 64) {
    throw std::invalid_argument("tensor Gram form supports q <= 64");
  }
  std::vector<Mask> form;
  form.reserve(label_rows.size());
  for (const auto left : label_rows) {
    Mask row = 0;
    for (std::size_t right = 0; right < label_rows.size(); ++right) {
      if (dot(left, label_rows[right])) row |= Mask{1} << right;
    }
    form.push_back(row);
  }
  return form;
}

std::uint32_t tensor_bilinear_rank(std::span<const Mask> label_rows) {
  return row_rank(tensor_bilinear_form_rows(label_rows));
}

bool tensor_bilinear_is_alternating(std::span<const Mask> label_rows) {
  return std::all_of(
      label_rows.begin(), label_rows.end(),
      [](Mask row) { return (std::popcount(row) & 1U) == 0; });
}

std::uint64_t nondegenerate_tensor_hyperplane_count(
    std::span<const Mask> label_rows) {
  const auto q = static_cast<std::uint32_t>(label_rows.size());
  if (q < 2 || q >= 63) {
    throw std::invalid_argument(
        "tensor hyperplane enumeration requires 2 <= q < 63");
  }
  std::uint64_t count = 0;
  const auto functional_count = std::uint64_t{1} << q;
  for (Mask functional = 1; functional < functional_count; ++functional) {
    const std::array<Mask, 1> constraint{functional};
    const auto kernel = nullspace_basis(constraint, q);
    std::vector<Mask> restricted_rows;
    restricted_rows.reserve(kernel.size());
    for (const auto coefficients : kernel) {
      restricted_rows.push_back(
          linear_combination(coefficients, label_rows));
    }
    if (tensor_radical_basis(restricted_rows).empty()) ++count;
  }
  return count;
}

bool is_primitive_exceptional_tensor(std::span<const Mask> label_rows) {
  const auto q = static_cast<std::uint32_t>(label_rows.size());
  if (q < 2) {
    throw std::invalid_argument("primitive tensor test requires q >= 2");
  }
  const auto bilinear_rank = tensor_bilinear_rank(label_rows);
  const bool exceptional_sector =
      bilinear_rank <= 1 ||
      ((q & 1U) == 0 && bilinear_rank == q &&
       tensor_bilinear_is_alternating(label_rows));
  const auto nondegenerate = q <= 7
                                 ? cubic_tensor_radical_basis_from_signature(
                                       q, cubic_tensor_word(label_rows))
                                       .empty()
                                 : tensor_radical_basis(label_rows).empty();
  return exceptional_sector && nondegenerate &&
         nondegenerate_tensor_hyperplane_count(label_rows) == 0;
}

LabelSpace rank_at_most_one_tensor_sector(const LabelSpace& label_space) {
  LabelSpace result = label_space;
  const auto gram = tensor_bilinear_form_rows(label_space.quotient_basis);
  Mask parity = 0;
  for (std::size_t index = 0; index < gram.size(); ++index) {
    if ((gram[index] >> index) & 1U) parity |= Mask{1} << index;
  }
  std::vector<Mask> alternating = gram;
  for (std::size_t index = 0; index < alternating.size(); ++index) {
    if ((parity >> index) & 1U) alternating[index] ^= parity;
    if ((alternating[index] >> index) & 1U) {
      throw std::logic_error("rank-one tensor-sector form is not alternating");
    }
  }
  result.bilinear_form_rows.push_back(std::move(alternating));
  return result;
}

LabelSpace even_weight_tensor_sector(const LabelSpace& label_space) {
  const auto width = label_space.quotient_dimension();
  Mask parity = 0;
  for (std::size_t index = 0; index < label_space.quotient_basis.size();
       ++index) {
    if (std::popcount(label_space.quotient_basis[index]) & 1U) {
      parity |= Mask{1} << index;
    }
  }
  if (parity == 0) return label_space;
  const std::array<Mask, 1> parity_constraint{parity};
  const auto kernel = nullspace_basis(parity_constraint, width);

  LabelSpace result;
  result.points = label_space.points;
  result.ambient_dimension = label_space.ambient_dimension;
  result.coordinate_masks = label_space.coordinate_masks;
  result.stabilizer_basis = label_space.stabilizer_basis;
  result.quotient_basis.reserve(kernel.size());
  for (const auto lift : kernel) {
    result.quotient_basis.push_back(
        linear_combination(lift, label_space.quotient_basis));
  }
  result.bilinear_form_rows.reserve(
      label_space.bilinear_form_rows.size());
  for (const auto& form : label_space.bilinear_form_rows) {
    std::vector<Mask> induced;
    induced.reserve(kernel.size());
    for (const auto left : kernel) {
      const auto contraction = linear_combination(left, form);
      Mask row = 0;
      for (std::size_t right = 0; right < kernel.size(); ++right) {
        if (dot(contraction, kernel[right])) row |= Mask{1} << right;
      }
      induced.push_back(row);
    }
    result.bilinear_form_rows.push_back(std::move(induced));
  }
  for (const auto row : result.quotient_basis) {
    if (std::popcount(row) & 1U) {
      throw std::logic_error("even-weight tensor sector contains an odd row");
    }
  }
  return result;
}

std::uint32_t intrinsic_tensor_dimension(std::span<const Mask> label_rows) {
  return static_cast<std::uint32_t>(label_rows.size()) -
         static_cast<std::uint32_t>(tensor_radical_basis(label_rows).size());
}

DistanceContext::DistanceContext(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension)
    : points_(points.begin(), points.end()),
      ambient_dimension_(ambient_dimension) {
  if (points_.size() > 64) {
    throw std::invalid_argument("distance context supports at most 64 columns");
  }
  if (ambient_dimension_ >= 32) {
    throw std::invalid_argument("distance context supports h < 32");
  }
  const std::size_t point_space = std::size_t{1} << ambient_dimension_;
  std::vector<std::int32_t> point_index(point_space, -1);
  for (std::size_t index = 0; index < points_.size(); ++index) {
    if (points_[index] >= point_space || point_index[points_[index]] != -1) {
      throw std::invalid_argument("distance support is invalid or repeated");
    }
    point_index[points_[index]] = static_cast<std::int32_t>(index);
  }
  for (std::size_t left = 0; left < points_.size(); ++left) {
    for (std::size_t right = left + 1; right < points_.size(); ++right) {
      const auto third = point_index[points_[left] ^ points_[right]];
      if (third <= static_cast<std::int32_t>(right)) continue;
      weight_three_subsets_.push_back(
          (Mask{1} << left) | (Mask{1} << right) | (Mask{1} << third));
    }
  }
}

std::uint64_t DistanceContext::weight_three_error_coefficient(
    std::span<const Mask> label_rows) const {
  std::uint64_t coefficient = 0;
  for (const Mask subset : weight_three_subsets_) {
    bool logical_nonzero = false;
    for (const Mask row : label_rows) {
      if (std::popcount(row & subset) & 1U) {
        logical_nonzero = true;
        break;
      }
    }
    coefficient += logical_nonzero;
  }
  return coefficient;
}

std::uint64_t DistanceContext::weight_four_error_coefficient(
    std::span<const Mask> label_rows) const {
  if (!weight_four_ready_) {
    std::unordered_map<std::uint32_t, std::vector<Mask>> pairs_by_xor;
    for (std::size_t left = 0; left < points_.size(); ++left) {
      for (std::size_t right = left + 1; right < points_.size(); ++right) {
        pairs_by_xor[points_[left] ^ points_[right]].push_back(
            (Mask{1} << left) | (Mask{1} << right));
      }
    }
    for (const auto& [syndrome, pairs] : pairs_by_xor) {
      (void)syndrome;
      for (std::size_t left = 0; left < pairs.size(); ++left) {
        for (std::size_t right = left + 1; right < pairs.size(); ++right) {
          if ((pairs[left] & pairs[right]) == 0) {
            weight_four_subsets_.push_back(pairs[left] | pairs[right]);
          }
        }
      }
    }
    std::sort(weight_four_subsets_.begin(), weight_four_subsets_.end());
    weight_four_subsets_.erase(
        std::unique(weight_four_subsets_.begin(), weight_four_subsets_.end()),
        weight_four_subsets_.end());
    weight_four_ready_ = true;
  }

  std::uint64_t coefficient = 0;
  for (const Mask subset : weight_four_subsets_) {
    bool logical_nonzero = false;
    for (const Mask row : label_rows) {
      if (std::popcount(row & subset) & 1U) {
        logical_nonzero = true;
        break;
      }
    }
    coefficient += logical_nonzero;
  }
  return coefficient;
}

DistanceResult DistanceContext::evaluate(
    std::span<const Mask> label_rows,
    std::uint32_t maximum_dp_dimension) const {
  const auto coefficient = weight_three_error_coefficient(label_rows);
  if (coefficient != 0) return DistanceResult{3U, coefficient};
  const auto weight_four_coefficient = weight_four_error_coefficient(label_rows);
  if (weight_four_coefficient != 0) {
    return DistanceResult{4U, weight_four_coefficient};
  }

  const auto q = static_cast<std::uint32_t>(label_rows.size());
  if (q >= 32 || ambient_dimension_ >= 32 || q + ambient_dimension_ >= 32) {
    throw std::invalid_argument("distance syndrome representation exceeds 31 bits");
  }
  const auto total_dimension = q + ambient_dimension_;
  if (total_dimension > maximum_dp_dimension) {
    throw std::invalid_argument("exact distance requires a larger DP limit");
  }

  std::vector<std::uint32_t> labels(points_.size(), 0);
  for (std::size_t index = 0; index < points_.size(); ++index) {
    for (std::uint32_t logical = 0; logical < q; ++logical) {
      labels[index] |= ((label_rows[logical] >> index) & 1U) << logical;
    }
  }
  const std::size_t syndrome_count = std::size_t{1} << total_dimension;
  const auto infinity = static_cast<std::uint32_t>(points_.size() + 1);
  std::vector<std::uint32_t> distances(syndrome_count, infinity);
  std::vector<std::uint64_t> counts(syndrome_count, 0);
  distances[0] = 0;
  counts[0] = 1;
  for (std::size_t index = 0; index < points_.size(); ++index) {
    const std::uint32_t column = labels[index] | (points_[index] << q);
    auto next_distances = distances;
    auto next_counts = counts;
    for (std::size_t syndrome = 0; syndrome < syndrome_count; ++syndrome) {
      if (distances[syndrome] == infinity) continue;
      const auto target = syndrome ^ column;
      const auto candidate = distances[syndrome] + 1;
      if (candidate < next_distances[target]) {
        next_distances[target] = candidate;
        next_counts[target] = counts[syndrome];
      } else if (candidate == next_distances[target]) {
        next_counts[target] += counts[syndrome];
      }
    }
    distances.swap(next_distances);
    counts.swap(next_counts);
  }
  std::uint32_t best = infinity;
  std::uint64_t best_count = 0;
  const std::size_t logical_count = std::size_t{1} << q;
  for (std::size_t target = 1; target < logical_count; ++target) {
    if (distances[target] < best) {
      best = distances[target];
      best_count = counts[target];
    } else if (distances[target] == best) {
      best_count += counts[target];
    }
  }
  if (best == infinity) return DistanceResult{std::nullopt, 0};
  return DistanceResult{best, best_count};
}

DistanceResult z_distance_and_error_coefficient(
    std::span<const std::uint32_t> points,
    std::uint32_t ambient_dimension,
    std::span<const Mask> label_rows,
    std::uint32_t maximum_dp_dimension) {
  return DistanceContext(points, ambient_dimension)
      .evaluate(label_rows, maximum_dp_dimension);
}

IsotropicSummary summarize_isotropic_subspaces(
    const LabelSpace& label_space, std::uint32_t dimension) {
  IsotropicSummary summary;
  summary.subspace_count = enumerate_totally_isotropic_subspaces(
      label_space, dimension,
      [&](std::span<const Mask> quotient_rows, std::span<const Mask> label_rows) {
        fnv_word(summary.tensor_checksum, dimension);
        for (const Mask row : quotient_rows) fnv_word(summary.tensor_checksum, row);
        for (const auto word : cubic_tensor_words(label_rows)) {
          fnv_word(summary.tensor_checksum, word);
        }
      });
  return summary;
}

std::string hex_mask(Mask value) {
  std::ostringstream output;
  output << std::hex << value;
  return output.str();
}

}  // namespace utsp
