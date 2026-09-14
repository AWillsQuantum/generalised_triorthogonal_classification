#pragma once

#include <cstdint>
#include <set>
#include <span>
#include <string_view>
#include <vector>

namespace utsp {

enum class PrimitiveGramSector {
  zero,
  rank_one,
  full_alternating,
};

[[nodiscard]] std::string_view primitive_gram_sector_name(
    PrimitiveGramSector sector);

struct PrimitiveTensorOrbit {
  std::uint64_t alternating_representative = 0;
  std::uint64_t standard_tensor_signature = 0;
  std::uint64_t orbit_size = 0;
  std::uint64_t stabilizer_order = 0;
  std::uint64_t nondegenerate_hyperplanes = 0;
};

struct TensorHyperplaneBucket {
  std::uint64_t nondegenerate_hyperplanes = 0;
  std::uint64_t orbit_count = 0;
  std::uint64_t tensor_count = 0;
};

struct PrimitiveTensorSectorCensus {
  std::uint32_t logical_qubits = 0;
  PrimitiveGramSector sector = PrimitiveGramSector::zero;
  std::uint64_t repeated_signature = 0;
  std::uint32_t alternating_dimension = 0;
  std::uint32_t action_generator_count = 0;
  std::uint64_t gram_stabilizer_order = 0;
  std::uint64_t tensor_count = 0;
  std::uint64_t tensor_orbit_count = 0;
  std::uint64_t nondegenerate_tensor_count = 0;
  std::uint64_t nondegenerate_tensor_orbit_count = 0;
  std::uint64_t primitive_tensor_count = 0;
  std::uint64_t primitive_tensor_orbit_count = 0;
  std::uint32_t generated_symplectic_transvection_labels = 0;
  std::vector<TensorHyperplaneBucket> nondegenerate_hyperplane_histogram;
  std::vector<PrimitiveTensorOrbit> primitive_orbits;
};

struct Q7AlternatingFormOrbit {
  std::uint32_t authority_id = 0;
  std::uint64_t alternating_signature = 0;
  std::uint64_t stabilizer_order = 0;
  std::uint64_t orbit_size = 0;
  std::uint32_t radical_dimension = 0;
  std::uint64_t nondegenerate_hyperplanes = 0;
};

struct Q7PrimitiveTensorOrbit {
  PrimitiveGramSector sector = PrimitiveGramSector::zero;
  std::uint32_t authority_id = 0;
  std::uint32_t rank_one_functional = 0;
  std::uint64_t standard_tensor_signature = 0;
  std::vector<std::uint64_t> direct_canonical_key;
  std::uint64_t represented_authority_pairs = 0;
};

struct Q7PrimitiveTensorAuthorityCensus {
  std::uint64_t general_linear_group_order = 0;
  std::uint64_t alternating_tensor_count = 0;
  std::uint64_t alternating_orbit_mass = 0;
  std::uint64_t rank_one_pairs_tested = 0;
  std::uint64_t rank_one_nondegenerate_pairs = 0;
  std::uint64_t rank_one_primitive_pairs = 0;
  std::vector<Q7AlternatingFormOrbit> alternating_orbits;
  std::vector<Q7PrimitiveTensorOrbit> primitive_orbits;
};

struct Q8AlternatingFormOrbit {
  std::uint32_t authority_id = 0;
  std::uint64_t alternating_state = 0;
  std::uint64_t stabilizer_order = 0;
  std::uint64_t orbit_size = 0;
  std::uint32_t radical_dimension = 0;
};

struct Q8PrimitiveTensorCandidate {
  PrimitiveGramSector sector = PrimitiveGramSector::zero;
  std::uint32_t authority_id = 0;
  std::uint32_t rank_one_functional = 0;
  std::uint64_t alternating_state = 0;
};

struct Q8PrimitiveTensorAuthorityScreen {
  std::uint64_t general_linear_group_order = 0;
  std::uint64_t alternating_tensor_count = 0;
  std::uint64_t alternating_orbit_mass = 0;
  std::uint64_t rank_one_pairs_tested = 0;
  std::uint64_t rank_one_nondegenerate_pairs = 0;
  std::uint64_t rank_one_primitive_pairs = 0;
  bool full_alternating_sector_theorem_excluded = false;
  std::vector<Q8AlternatingFormOrbit> alternating_orbits;
  std::vector<Q8PrimitiveTensorCandidate> primitive_candidates;
};

struct TensorRestrictionProfile {
  std::uint32_t target_dimension = 0;
  bool covers_all_primitive_orbits = false;
  bool uses_q7_primitive_family_recognizer = false;
  std::vector<std::set<std::vector<std::uint64_t>>> keys_by_dimension;

  [[nodiscard]] bool allows(
      std::uint32_t dimension,
      std::span<const std::uint64_t> canonical_key) const;

  [[nodiscard]] bool allows_standard_signature(
      std::uint32_t dimension,
      std::uint64_t standard_tensor_signature) const;
};

[[nodiscard]] std::vector<PrimitiveGramSector> primitive_gram_sectors(
    std::uint32_t logical_qubits);

// Enumerate the alternating-part orbits under the exact stabilizer of a
// canonical exceptional Gram form.  Combined with the primitive-sector
// theorem, these are all intrinsic tensors not inherited from a
// nondegenerate hyperplane.
[[nodiscard]] PrimitiveTensorSectorCensus census_primitive_tensor_sector(
    std::uint32_t logical_qubits,
    PrimitiveGramSector sector);

// The complete Cohen--Helminck dimension-seven alternating-form authority,
// independently audited by the published stabilizer orders.  Rank-one Gram
// tensors are the unique pairs T=A+p^3, so testing the 12 authority forms and
// all 127 nonzero p covers every GL(7,2)-orbit without enumerating 2^35 forms.
[[nodiscard]] Q7PrimitiveTensorAuthorityCensus
census_q7_primitive_tensors_from_authority(std::uint32_t workers = 1);

// Complete zero/rank-one q=8 screen over the 32 Hora--Pudlak authority
// forms.  Rank-one entries are candidate authority pairs until stabilizer
// orbits of their functionals have been certified.
[[nodiscard]] Q8PrimitiveTensorAuthorityScreen
screen_q8_primitive_tensors_from_authority();

// Exact GL-invariant keys of every nonzero-dimensional restriction of the
// supplied target tensors.  A partial support subspace can extend to a
// target only if its induced tensor key occurs at its current dimension.
[[nodiscard]] TensorRestrictionProfile build_tensor_restriction_profile(
    std::uint32_t logical_qubits,
    std::span<const std::uint64_t> standard_tensor_signatures,
    std::uint32_t minimum_dimension = 1,
    std::uint32_t workers = 1);

[[nodiscard]] bool has_nondegenerate_tensor_restriction(
    std::uint32_t logical_qubits,
    std::uint64_t standard_tensor_signature,
    std::uint32_t restriction_dimension);

}  // namespace utsp
