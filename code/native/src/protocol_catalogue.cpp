#include "utsp/protocol_catalogue.hpp"

#include "utsp/known_outputs.hpp"
#include "utsp/marked_code_canonical.hpp"
#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <exception>
#include <map>
#include <set>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace utsp {
namespace {

struct TensorWordsHash {
  [[nodiscard]] std::size_t operator()(
      const std::vector<std::uint64_t>& words) const noexcept {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for (const auto word : words) {
      hash ^= word;
      hash *= 0x100000001b3ULL;
      hash ^= word >> 32U;
      hash *= 0x100000001b3ULL;
    }
    return static_cast<std::size_t>(hash);
  }
};

[[nodiscard]] const KnownOutputClassifier& known_output_classifier() {
  // Construction builds all q<=4 completion and GL-orbit tables.  The tables
  // are immutable, so one process-wide instance is both safe and substantially
  // cheaper than rebuilding them in every support microchunk.
  static const KnownOutputClassifier classifier;
  return classifier;
}

struct TensorAnalysis {
  std::uint32_t radical_dimension = 0;
  OutputTensorKey canonical_key;
  std::vector<std::uint32_t> canonical_basis;
  std::string known_class_id;
};

struct FrontierPartition {
  OutputTensorKey output;
  std::uint32_t distance = 0;

  [[nodiscard]] bool operator<(const FrontierPartition& other) const {
    return std::tie(output, distance) <
           std::tie(other.output, other.distance);
  }
};

using TensorCache = std::unordered_map<
    std::vector<std::uint64_t>, TensorAnalysis, TensorWordsHash>;
using FrontierMap =
    std::map<FrontierPartition, std::vector<ProtocolFrontierWitness>>;
using OutputOrbitSet = std::set<OutputTensorKey>;
using RadicalDimensionCounts = std::map<std::uint32_t, std::uint64_t>;

[[nodiscard]] bool witness_tie_less(
    const ProtocolFrontierWitness& left,
    const ProtocolFrontierWitness& right) {
  return std::tie(
             left.source_id,
             left.source_family,
             left.source_space_length,
             left.parity_case,
             left.origin,
             left.points,
             left.canonical_logical_rows,
             left.record_index) <
         std::tie(
             right.source_id,
             right.source_family,
             right.source_space_length,
             right.parity_case,
             right.origin,
             right.points,
             right.canonical_logical_rows,
             right.record_index);
}

[[nodiscard]] std::vector<Mask> transform_label_rows(
    std::span<const Mask> label_rows,
    std::span<const std::uint32_t> basis);

void update_frontier(
    std::vector<ProtocolFrontierWitness>& frontier,
    const TensorAnalysis& tensor,
    const DistanceResult& distance,
    const ProtocolSupport& support,
    std::span<const Mask> label_rows,
    std::uint32_t logical_qubits) {
  const auto protocol_length =
      static_cast<std::uint32_t>(support.points.size());
  const auto space_footprint = support.ambient_dimension + logical_qubits;
  for (auto& existing : frontier) {
    if (existing.protocol_length == protocol_length &&
        existing.space_footprint == space_footprint) {
      ProtocolFrontierWitness candidate;
      candidate.output = tensor.canonical_key;
      candidate.known_class_id = tensor.known_class_id;
      candidate.distance = *distance.distance;
      candidate.error_coefficient = distance.error_coefficient;
      candidate.protocol_length = protocol_length;
      candidate.space_footprint = space_footprint;
      candidate.record_index = support.record_index;
      candidate.source_index = support.source_index;
      candidate.source_id = support.source_id;
      candidate.source_family = support.source_family;
      candidate.source_space_length = support.source_space_length;
      candidate.ambient_dimension = support.ambient_dimension;
      candidate.parity_case = support.parity_case;
      candidate.origin = support.origin;
      candidate.points = support.points;
      candidate.canonical_logical_rows = transform_label_rows(
          label_rows, tensor.canonical_basis);
      if (witness_tie_less(candidate, existing)) existing = std::move(candidate);
      return;
    }
    if (existing.protocol_length <= protocol_length &&
        existing.space_footprint <= space_footprint) {
      return;
    }
  }
  frontier.erase(
      std::remove_if(
          frontier.begin(), frontier.end(),
          [&](const ProtocolFrontierWitness& existing) {
            return protocol_length <= existing.protocol_length &&
                   space_footprint <= existing.space_footprint;
          }),
      frontier.end());
  ProtocolFrontierWitness candidate;
  candidate.output = tensor.canonical_key;
  candidate.known_class_id = tensor.known_class_id;
  candidate.distance = *distance.distance;
  candidate.error_coefficient = distance.error_coefficient;
  candidate.protocol_length = protocol_length;
  candidate.space_footprint = space_footprint;
  candidate.record_index = support.record_index;
  candidate.source_index = support.source_index;
  candidate.source_id = support.source_id;
  candidate.source_family = support.source_family;
  candidate.source_space_length = support.source_space_length;
  candidate.ambient_dimension = support.ambient_dimension;
  candidate.parity_case = support.parity_case;
  candidate.origin = support.origin;
  candidate.points = support.points;
  candidate.canonical_logical_rows = transform_label_rows(
      label_rows, tensor.canonical_basis);
  frontier.push_back(std::move(candidate));
}

[[nodiscard]] std::vector<Mask> transform_label_rows(
    std::span<const Mask> label_rows,
    std::span<const std::uint32_t> basis) {
  std::vector<Mask> transformed;
  transformed.reserve(basis.size());
  for (const auto coefficients : basis) {
    transformed.push_back(linear_combination(coefficients, label_rows));
  }
  return transformed;
}

void update_witness_frontier(
    std::vector<ProtocolFrontierWitness>& frontier,
    ProtocolFrontierWitness candidate) {
  for (auto& existing : frontier) {
    if (existing.protocol_length == candidate.protocol_length &&
        existing.space_footprint == candidate.space_footprint) {
      if (witness_tie_less(candidate, existing)) {
        existing = std::move(candidate);
      }
      return;
    }
    if (existing.protocol_length <= candidate.protocol_length &&
        existing.space_footprint <= candidate.space_footprint) {
      return;
    }
  }
  frontier.erase(
      std::remove_if(
          frontier.begin(),
          frontier.end(),
          [&](const ProtocolFrontierWitness& existing) {
            return candidate.protocol_length <= existing.protocol_length &&
                   candidate.space_footprint <= existing.space_footprint;
          }),
      frontier.end());
  frontier.push_back(std::move(candidate));
}

}  // namespace

bool OutputTensorKey::operator<(const OutputTensorKey& other) const {
  return std::tie(logical_qubits, words) <
         std::tie(other.logical_qubits, other.words);
}

std::size_t QuotientBasisCacheView::size() const {
  return offsets.empty() ? 0 : offsets.size() - 1;
}

std::span<const Mask> QuotientBasisCacheView::basis(
    std::size_t index) const {
  if (!provided() || index >= size()) {
    throw std::out_of_range("cached quotient-basis index is out of range");
  }
  const auto begin = offsets[index];
  const auto end = offsets[index + 1];
  if (begin > end || end > rows.size()) {
    throw std::invalid_argument("cached quotient-basis offsets are invalid");
  }
  return rows.subspan(
      static_cast<std::size_t>(begin), static_cast<std::size_t>(end - begin));
}

QuotientBasisCacheView QuotientBasisCacheView::subview(
    std::size_t begin, std::size_t count) const {
  if (!provided()) return {};
  if (begin > size() || count > size() - begin) {
    throw std::out_of_range("cached quotient-basis subview is out of range");
  }
  return QuotientBasisCacheView{offsets.subspan(begin, count + 1), rows};
}

ProtocolCatalogueResult classify_protocol_frontier(
    std::span<const ProtocolSupport> supports,
    std::uint32_t logical_qubits,
    std::uint32_t minimum_distance,
    std::uint32_t maximum_distance_dp_dimension,
    ProtocolEnumerationMode enumeration_mode,
    std::uint32_t raw_maximum_quotient_dimension,
    std::uint32_t orbit_workers,
    std::uint32_t support_workers,
    std::uint32_t zero_isotropic_workers,
    std::uint32_t minimum_quotient_dimension,
    std::uint32_t maximum_quotient_dimension,
    bool collect_positive_supports,
    QuotientBasisCacheView cached_quotient_bases) {
  if (logical_qubits < 1 || logical_qubits > 8) {
    throw std::invalid_argument(
        "direct protocol classification currently supports q=1,...,8");
  }
  if ((enumeration_mode ==
           ProtocolEnumerationMode::q3_predecessor_target_marked_code_orbits ||
       enumeration_mode ==
           ProtocolEnumerationMode::q3_predecessor_target_existence) &&
      logical_qubits != 3) {
    throw std::invalid_argument(
        "q3 predecessor target enumeration requires q=3");
  }
  if (enumeration_mode ==
          ProtocolEnumerationMode::q7_primitive_marked_code_orbits &&
      logical_qubits != 7) {
    throw std::invalid_argument(
        "q7 primitive marked-code enumeration requires q=7");
  }
  if (enumeration_mode ==
          ProtocolEnumerationMode::q5_q3_chain_cover_marked_code_orbits &&
      logical_qubits != 5) {
    throw std::invalid_argument(
        "q5 q3-chain-cover marked-code enumeration requires q=5");
  }
  if (enumeration_mode ==
          ProtocolEnumerationMode::q5_q4_hitting_set_marked_code_orbits &&
      logical_qubits != 5) {
    throw std::invalid_argument(
        "q5 q4-hitting-set marked-code enumeration requires q=5");
  }
  if (orbit_workers == 0 || support_workers == 0) {
    throw std::invalid_argument("protocol worker counts must be positive");
  }
  if (zero_isotropic_workers != 0 &&
      enumeration_mode != ProtocolEnumerationMode::raw_subspaces) {
    throw std::invalid_argument(
        "zero-isotropic proof requires raw subspace enumeration");
  }
  if (minimum_quotient_dimension > maximum_quotient_dimension) {
    throw std::invalid_argument("invalid quotient-dimension interval");
  }
  if (cached_quotient_bases.provided() &&
      cached_quotient_bases.size() != supports.size()) {
    throw std::invalid_argument(
        "cached quotient-basis count does not match support count");
  }
  if (support_workers > 1 && supports.size() > 1) {
    const auto worker_count = static_cast<std::size_t>(
        std::min<std::uint64_t>(support_workers, supports.size()));
    // Higher-q marked-code supports have a heavy-tailed cost distribution.
    // Target small dynamic chunks so a few rare high-quotient supports do not
    // serialize an otherwise parallel block, while capping scheduler state.
    const auto task_count = [&] {
      if (logical_qubits < 3) {
        return std::min(supports.size(), worker_count * std::size_t{8});
      }
      constexpr std::size_t target_supports_per_task = 16;
      constexpr std::size_t maximum_tasks_per_worker = 16384;
      const auto target_tasks =
          (supports.size() + target_supports_per_task - 1) /
          target_supports_per_task;
      return std::min(
          supports.size(),
          std::max(
              worker_count,
              std::min(
                  target_tasks, worker_count * maximum_tasks_per_worker)));
    }();
    std::vector<ProtocolCatalogueResult> partial(task_count);
    std::vector<std::exception_ptr> errors(task_count);
    std::vector<std::size_t> task_order(task_count);
    for (std::size_t task = 0; task < task_count; ++task) {
      task_order[task] = task;
    }
    if (logical_qubits >= 3 && cached_quotient_bases.provided()) {
      struct TaskCost {
        std::uint32_t maximum_dimension = 0;
        std::uint64_t exponential_weight = 0;
      };
      std::vector<TaskCost> task_costs(task_count);
      for (std::size_t task = 0; task < task_count; ++task) {
        const auto begin = supports.size() * task / task_count;
        const auto end = supports.size() * (task + 1) / task_count;
        auto& cost = task_costs[task];
        for (auto index = begin; index < end; ++index) {
          const auto dimension = static_cast<std::uint32_t>(
              cached_quotient_bases.basis(index).size());
          cost.maximum_dimension =
              std::max(cost.maximum_dimension, dimension);
          cost.exponential_weight +=
              std::uint64_t{1} << std::min(dimension, std::uint32_t{55});
        }
      }
      std::stable_sort(
          task_order.begin(),
          task_order.end(),
          [&](const auto left, const auto right) {
            const auto& left_cost = task_costs[left];
            const auto& right_cost = task_costs[right];
            return std::tie(
                       left_cost.maximum_dimension,
                       left_cost.exponential_weight) >
                   std::tie(
                       right_cost.maximum_dimension,
                       right_cost.exponential_weight);
          });
    }
    std::atomic<std::size_t> next_task{0};
    std::vector<std::thread> threads;
    threads.reserve(worker_count);
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
      threads.emplace_back([&] {
        while (true) {
          const auto scheduled_task = next_task.fetch_add(1);
          if (scheduled_task >= task_count) return;
          const auto task = task_order[scheduled_task];
          const auto begin = supports.size() * task / task_count;
          const auto end = supports.size() * (task + 1) / task_count;
          try {
            partial[task] = classify_protocol_frontier(
                supports.subspan(begin, end - begin),
                logical_qubits,
                minimum_distance,
                maximum_distance_dp_dimension,
                enumeration_mode,
                raw_maximum_quotient_dimension,
                orbit_workers,
                1,
                zero_isotropic_workers,
                minimum_quotient_dimension,
                maximum_quotient_dimension,
                collect_positive_supports,
                cached_quotient_bases.subview(begin, end - begin));
          } catch (...) {
            errors[task] = std::current_exception();
            next_task.store(task_count);
            return;
          }
        }
      });
    }
    for (auto& thread : threads) thread.join();
    for (const auto& error : errors) {
      if (error) std::rethrow_exception(error);
    }

    ProtocolCatalogueResult merged;
    merged.logical_qubits = logical_qubits;
    merged.minimum_distance = minimum_distance;
    merged.enumeration_mode = enumeration_mode;
    merged.raw_maximum_quotient_dimension =
        raw_maximum_quotient_dimension;
    merged.zero_isotropic_proof = zero_isotropic_workers != 0;
    merged.zero_isotropic_workers = zero_isotropic_workers;
    std::set<std::vector<std::uint64_t>> raw_tensor_keys;
    std::set<OutputTensorKey> output_keys;
    std::map<FrontierPartition, std::vector<ProtocolFrontierWitness>>
        frontiers;
    for (std::size_t task = 0; task < partial.size(); ++task) {
      auto& child = partial[task];
      const auto begin = supports.size() * task / task_count;
      merged.isotropic_subspace_statistics_exact &=
          child.isotropic_subspace_statistics_exact;
      merged.radical_statistics_exact &= child.radical_statistics_exact;
      merged.supports_processed += child.supports_processed;
      merged.quotient_dimension_filtered_supports +=
          child.quotient_dimension_filtered_supports;
      merged.eligible_supports += child.eligible_supports;
      merged.raw_enumerated_supports += child.raw_enumerated_supports;
      merged.marked_orbit_enumerated_supports +=
          child.marked_orbit_enumerated_supports;
      merged.isotropic_subspaces += child.isotropic_subspaces;
      merged.nondegenerate_subspaces += child.nondegenerate_subspaces;
      merged.marked_code_orbits += child.marked_code_orbits;
      merged.exact_distance_calls += child.exact_distance_calls;
      for (const auto& [dimension, count] :
           child.radical_dimension_counts) {
        merged.radical_dimension_counts[dimension] += count;
      }
      raw_tensor_keys.insert(
          child.raw_tensor_keys.begin(), child.raw_tensor_keys.end());
      output_keys.insert(
          child.canonical_output_keys.begin(),
          child.canonical_output_keys.end());
      for (auto& positive : child.positive_supports) {
        positive.support_offset += begin;
        merged.positive_supports.push_back(std::move(positive));
      }
      for (auto& witness : child.pareto_witnesses) {
        const FrontierPartition partition{witness.output, witness.distance};
        update_witness_frontier(
            frontiers[partition],
            std::move(witness));
      }
    }
    merged.raw_tensor_keys.assign(
        raw_tensor_keys.begin(), raw_tensor_keys.end());
    merged.raw_tensor_forms = merged.raw_tensor_keys.size();
    merged.canonical_output_keys.assign(output_keys.begin(), output_keys.end());
    merged.canonical_output_orbits = merged.canonical_output_keys.size();
    for (auto& [partition, witnesses] : frontiers) {
      (void)partition;
      std::sort(
          witnesses.begin(),
          witnesses.end(),
          [](const auto& left, const auto& right) {
            return std::tie(left.protocol_length, left.space_footprint) <
                   std::tie(right.protocol_length, right.space_footprint);
          });
      for (auto& witness : witnesses) {
        merged.pareto_witnesses.push_back(std::move(witness));
      }
    }
    return merged;
  }
  ProtocolCatalogueResult result;
  result.logical_qubits = logical_qubits;
  result.minimum_distance = minimum_distance;
  result.enumeration_mode = enumeration_mode;
  result.raw_maximum_quotient_dimension =
      raw_maximum_quotient_dimension;
  result.zero_isotropic_proof = zero_isotropic_workers != 0;
  result.zero_isotropic_workers = zero_isotropic_workers;
  result.supports_processed = supports.size();
  TensorCache tensor_cache;
  FrontierMap frontiers;
  OutputOrbitSet output_orbits;

  for (std::size_t support_offset = 0; support_offset < supports.size();
       ++support_offset) {
    const auto& support = supports[support_offset];
    std::span<const Mask> cached_basis;
    if (cached_quotient_bases.provided()) {
      cached_basis = cached_quotient_bases.basis(support_offset);
      const auto quotient_dimension =
          static_cast<std::uint32_t>(cached_basis.size());
      if (quotient_dimension < minimum_quotient_dimension ||
          quotient_dimension > maximum_quotient_dimension) {
        ++result.quotient_dimension_filtered_supports;
        continue;
      }
      if (quotient_dimension < logical_qubits) continue;
    }
    const auto label_space = !cached_quotient_bases.provided()
                                 ? build_logical_label_space(
                                       support.points,
                                       support.ambient_dimension,
                                       false,
                                       minimum_distance)
                                 : build_logical_label_space_from_quotient_basis(
                                       support.points,
                                       support.ambient_dimension,
                                       cached_basis,
                                       false);
    if (!cached_quotient_bases.provided() &&
        (label_space.quotient_dimension() < minimum_quotient_dimension ||
         label_space.quotient_dimension() > maximum_quotient_dimension)) {
      ++result.quotient_dimension_filtered_supports;
      continue;
    }
    if (!cached_quotient_bases.provided() &&
        label_space.quotient_dimension() < logical_qubits) {
      continue;
    }
    ++result.eligible_supports;
    const bool use_marked_orbits =
        enumeration_mode == ProtocolEnumerationMode::marked_code_orbits ||
        enumeration_mode == ProtocolEnumerationMode::
                                final_dimension_nondegenerate_marked_code_orbits ||
        enumeration_mode ==
            ProtocolEnumerationMode::q3_predecessor_target_marked_code_orbits ||
        enumeration_mode ==
            ProtocolEnumerationMode::q3_predecessor_target_existence ||
        enumeration_mode == ProtocolEnumerationMode::
                                q5_q3_chain_cover_marked_code_orbits ||
        enumeration_mode == ProtocolEnumerationMode::
                                q5_q4_hitting_set_marked_code_orbits ||
        enumeration_mode ==
            ProtocolEnumerationMode::q7_primitive_marked_code_orbits ||
        (enumeration_mode == ProtocolEnumerationMode::hybrid &&
         label_space.quotient_dimension() >
             raw_maximum_quotient_dimension);
    if (use_marked_orbits) {
      ++result.marked_orbit_enumerated_supports;
      result.isotropic_subspace_statistics_exact = false;
      result.radical_statistics_exact = false;
    } else {
      ++result.raw_enumerated_supports;
    }
    if (zero_isotropic_workers != 0) {
      const auto isotropic_count = count_totally_isotropic_subspaces_parallel(
          label_space, logical_qubits, zero_isotropic_workers);
      result.isotropic_subspaces += isotropic_count;
      if (isotropic_count != 0) {
        throw std::logic_error(
            "zero-isotropic proof encountered an isotropic subspace");
      }
      continue;
    }
    if (enumeration_mode ==
        ProtocolEnumerationMode::q3_predecessor_target_existence) {
      const auto existence =
          find_q3_predecessor_target_subspace(label_space, orbit_workers);
      result.marked_code_orbits +=
          existence.cs_seed_orbits + existence.primitive_target_orbits;
      if (existence.found) {
        std::vector<Mask> label_rows;
        label_rows.reserve(existence.witness_quotient_rows.size());
        for (const auto row : existence.witness_quotient_rows) {
          label_rows.push_back(
              linear_combination(row, label_space.quotient_basis));
        }
        const auto canonical =
            canonicalize_cubic_tensor_from_labels(label_rows);
        const OutputTensorKey output{
            logical_qubits, canonical.canonical_key_words};
        ++result.nondegenerate_subspaces;
        output_orbits.insert(output);
        if (collect_positive_supports) {
          result.positive_supports.push_back(ProtocolPositiveSupport{
              support_offset, support.record_index, {output}});
        }
      }
      continue;
    }
    const DistanceContext distance_context(
        support.points, support.ambient_dimension);
    OutputOrbitSet support_output_orbits;
    const auto process_subspace_into =
        [&](std::span<const Mask> quotient_rows,
            std::span<const Mask> label_rows,
            std::uint64_t orbit_size,
            const DistanceContext& local_distance_context,
            TensorCache& local_tensor_cache,
            FrontierMap& local_frontiers,
            OutputOrbitSet& local_output_orbits,
            OutputOrbitSet& local_support_output_orbits,
            std::uint64_t& local_nondegenerate_subspaces,
            std::uint64_t& local_exact_distance_calls,
            RadicalDimensionCounts& local_radical_dimension_counts) {
          (void)quotient_rows;
          const auto raw_tensor = cubic_tensor_words(label_rows);
          auto analysis = local_tensor_cache.find(raw_tensor);
          if (analysis == local_tensor_cache.end()) {
            TensorAnalysis computed;
            computed.radical_dimension = static_cast<std::uint32_t>(
                tensor_radical_basis(label_rows).size());
            if (computed.radical_dimension == 0) {
              const auto canonical =
                  canonicalize_cubic_tensor_from_labels(label_rows);
              computed.canonical_key = OutputTensorKey{
                  logical_qubits, canonical.canonical_key_words};
              computed.canonical_basis =
                  canonical.new_basis_in_old_coordinates;
              if (logical_qubits <= 4) {
                const auto* known = known_output_classifier().classify(
                    logical_qubits, cubic_tensor_word(label_rows));
                if (known == nullptr ||
                    known->intrinsic_logical_qubits != logical_qubits) {
                  throw std::logic_error(
                      "known-output table rejected a nondegenerate tensor");
                }
                computed.known_class_id = std::string(known->class_id);
              }
            }
            analysis = local_tensor_cache
                           .emplace(raw_tensor, std::move(computed))
                           .first;
          }
          const auto& tensor = analysis->second;
          if (tensor.radical_dimension != 0) {
            if (use_marked_orbits) {
              throw std::logic_error(
                  "complete marked-code enumeration returned a radical tensor");
            }
            local_radical_dimension_counts[tensor.radical_dimension] +=
                orbit_size;
            return;
          }

          local_nondegenerate_subspaces += orbit_size;
          ++local_exact_distance_calls;
          local_output_orbits.insert(tensor.canonical_key);
          if (collect_positive_supports) {
            local_support_output_orbits.insert(tensor.canonical_key);
          }
          const auto distance = local_distance_context.evaluate(
              label_rows, maximum_distance_dp_dimension);
          if (!distance.distance || *distance.distance < minimum_distance) {
            throw std::logic_error(
                "distance-filtered logical subspace failed its distance floor");
          }
          update_frontier(
              local_frontiers[FrontierPartition{
                  tensor.canonical_key, *distance.distance}],
              tensor,
              distance,
              support,
              label_rows,
              logical_qubits);
        };
    const auto process_subspace =
        [&](std::span<const Mask> quotient_rows,
            std::span<const Mask> label_rows,
            std::uint64_t orbit_size) {
          process_subspace_into(
              quotient_rows,
              label_rows,
              orbit_size,
              distance_context,
              tensor_cache,
              frontiers,
              output_orbits,
              support_output_orbits,
              result.nondegenerate_subspaces,
              result.exact_distance_calls,
              result.radical_dimension_counts);
        };

    if (use_marked_orbits) {
      const auto level = [&] {
        if (enumeration_mode ==
            ProtocolEnumerationMode::
                final_dimension_nondegenerate_marked_code_orbits) {
          return enumerate_dimension_filtered_nondegenerate_subspace_orbits(
              label_space,
              logical_qubits,
              logical_qubits,
              orbit_workers);
        }
        if (enumeration_mode ==
            ProtocolEnumerationMode::q3_predecessor_target_marked_code_orbits) {
          return enumerate_q3_predecessor_target_subspace_orbits(
              label_space, orbit_workers);
        }
        if (enumeration_mode == ProtocolEnumerationMode::
                                    q5_q3_chain_cover_marked_code_orbits) {
          return enumerate_q5_q3_chain_cover_subspace_orbits(
              label_space, orbit_workers);
        }
        if (enumeration_mode == ProtocolEnumerationMode::
                                    q5_q4_hitting_set_marked_code_orbits) {
          return enumerate_q5_q4_hitting_set_subspace_orbits(
              label_space, orbit_workers);
        }
        if (enumeration_mode ==
            ProtocolEnumerationMode::q7_primitive_marked_code_orbits) {
          constexpr std::array<std::uint64_t, 2> primitive_signatures{
              0x42010000000ULL,
              0x42010000001ULL,
          };
          return enumerate_target_tensor_subspace_orbits(
              label_space,
              logical_qubits,
              logical_qubits,
              primitive_signatures,
              3,
              orbit_workers);
        }
        return enumerate_complete_nondegenerate_subspace_orbits(
            label_space, logical_qubits, orbit_workers);
      }();
      if (level.representatives.size() != level.orbit_sizes.size()) {
        throw std::logic_error(
            "marked-code representatives and orbit sizes are misaligned");
      }
      result.marked_code_orbits += level.representatives.size();
      const auto evaluate_marked_orbit =
          [&](std::size_t index,
              const DistanceContext& local_distance_context,
              TensorCache& local_tensor_cache,
              FrontierMap& local_frontiers,
              OutputOrbitSet& local_output_orbits,
              OutputOrbitSet& local_support_output_orbits,
              std::uint64_t& local_nondegenerate_subspaces,
              std::uint64_t& local_exact_distance_calls,
              RadicalDimensionCounts& local_radical_dimension_counts) {
        const auto& quotient_rows = level.representatives[index];
        std::vector<Mask> label_rows;
        label_rows.reserve(quotient_rows.size());
        for (const auto row : quotient_rows) {
          label_rows.push_back(
              linear_combination(row, label_space.quotient_basis));
        }
        process_subspace_into(
            std::span<const Mask>(quotient_rows),
            std::span<const Mask>(label_rows),
            level.orbit_sizes[index],
            local_distance_context,
            local_tensor_cache,
            local_frontiers,
            local_output_orbits,
            local_support_output_orbits,
            local_nondegenerate_subspaces,
            local_exact_distance_calls,
            local_radical_dimension_counts);
      };
      const auto evaluation_workers = std::min<std::size_t>(
          orbit_workers, level.representatives.size());
      if (evaluation_workers <= 1) {
        for (std::size_t index = 0; index < level.representatives.size();
             ++index) {
          evaluate_marked_orbit(
              index,
              distance_context,
              tensor_cache,
              frontiers,
              output_orbits,
              support_output_orbits,
              result.nondegenerate_subspaces,
              result.exact_distance_calls,
              result.radical_dimension_counts);
        }
      } else {
        struct WorkerEvaluation {
          TensorCache tensor_cache;
          FrontierMap frontiers;
          OutputOrbitSet output_orbits;
          OutputOrbitSet support_output_orbits;
          RadicalDimensionCounts radical_dimension_counts;
          std::uint64_t nondegenerate_subspaces = 0;
          std::uint64_t exact_distance_calls = 0;
          std::exception_ptr error;
        };
        std::vector<WorkerEvaluation> evaluations(evaluation_workers);
        std::atomic<std::size_t> next_orbit{0};
        std::vector<std::thread> threads;
        threads.reserve(evaluation_workers);
        for (std::size_t worker = 0; worker < evaluation_workers; ++worker) {
          threads.emplace_back([&, worker] {
            DistanceContext local_distance_context(
                support.points, support.ambient_dimension);
            auto& evaluation = evaluations[worker];
            while (true) {
              const auto index = next_orbit.fetch_add(1);
              if (index >= level.representatives.size()) return;
              try {
                evaluate_marked_orbit(
                    index,
                    local_distance_context,
                    evaluation.tensor_cache,
                    evaluation.frontiers,
                    evaluation.output_orbits,
                    evaluation.support_output_orbits,
                    evaluation.nondegenerate_subspaces,
                    evaluation.exact_distance_calls,
                    evaluation.radical_dimension_counts);
              } catch (...) {
                evaluation.error = std::current_exception();
                next_orbit.store(level.representatives.size());
                return;
              }
            }
          });
        }
        for (auto& thread : threads) thread.join();
        for (const auto& evaluation : evaluations) {
          if (evaluation.error) std::rethrow_exception(evaluation.error);
        }
        for (auto& evaluation : evaluations) {
          result.nondegenerate_subspaces +=
              evaluation.nondegenerate_subspaces;
          result.exact_distance_calls += evaluation.exact_distance_calls;
          for (const auto& [dimension, count] :
               evaluation.radical_dimension_counts) {
            result.radical_dimension_counts[dimension] += count;
          }
          for (auto& [key, analysis] : evaluation.tensor_cache) {
            tensor_cache.try_emplace(key, std::move(analysis));
          }
          output_orbits.insert(
              evaluation.output_orbits.begin(),
              evaluation.output_orbits.end());
          support_output_orbits.insert(
              evaluation.support_output_orbits.begin(),
              evaluation.support_output_orbits.end());
          for (auto& [partition, witnesses] : evaluation.frontiers) {
            for (auto& witness : witnesses) {
              update_witness_frontier(
                  frontiers[partition], std::move(witness));
            }
          }
        }
      }
    } else {
      result.isotropic_subspaces += enumerate_totally_isotropic_subspaces(
          label_space,
          logical_qubits,
          [&](std::span<const Mask> quotient_rows,
              std::span<const Mask> label_rows) {
            process_subspace(quotient_rows, label_rows, 1);
           });
    }
    if (!support_output_orbits.empty()) {
      result.positive_supports.push_back(ProtocolPositiveSupport{
          support_offset,
          support.record_index,
          {support_output_orbits.begin(), support_output_orbits.end()},
      });
    }
  }

  result.raw_tensor_keys.reserve(tensor_cache.size());
  for (const auto& [key, analysis] : tensor_cache) {
    (void)analysis;
    result.raw_tensor_keys.push_back(key);
  }
  std::sort(result.raw_tensor_keys.begin(), result.raw_tensor_keys.end());
  result.raw_tensor_forms = result.raw_tensor_keys.size();
  result.canonical_output_keys.assign(
      output_orbits.begin(), output_orbits.end());
  result.canonical_output_orbits = result.canonical_output_keys.size();
  for (auto& [partition, witnesses] : frontiers) {
    (void)partition;
    std::sort(
        witnesses.begin(), witnesses.end(),
        [](const auto& left, const auto& right) {
          return std::tie(left.protocol_length, left.space_footprint) <
                 std::tie(right.protocol_length, right.space_footprint);
        });
    for (auto& witness : witnesses) {
      result.pareto_witnesses.push_back(std::move(witness));
    }
  }
  return result;
}

}  // namespace utsp
