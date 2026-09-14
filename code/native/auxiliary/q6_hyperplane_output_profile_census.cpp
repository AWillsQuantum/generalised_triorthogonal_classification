#include "utsp/tensor_canonical.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <compare>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

constexpr std::uint8_t kUnassigned = 0xff;

struct PackedAction25 {
  std::array<std::array<std::uint32_t, 32>, 5> chunks{};

  [[nodiscard]] std::uint32_t operator()(std::uint32_t signature) const {
    std::uint32_t result = 0;
    for (std::uint32_t chunk = 0; chunk < chunks.size(); ++chunk) {
      result ^= chunks[chunk][(signature >> (5 * chunk)) & 31U];
    }
    return result;
  }
};

struct RestrictionMap {
  std::array<std::vector<std::uint32_t>, 3> chunks;

  [[nodiscard]] std::uint32_t operator()(std::uint64_t signature) const {
    return chunks[0][signature & 0x3fffU] ^
           chunks[1][(signature >> 14) & 0x3fffU] ^
           chunks[2][(signature >> 28) & 0x1fffU];
  }
};

struct Q5Orbit {
  std::uint32_t canonical_key = 0;
  std::uint32_t representative = 0;
  std::uint64_t orbit_size = 0;
  std::uint32_t radical_dimension = 0;
  std::uint8_t nondegenerate_bit = kUnassigned;
};

struct Q5Census {
  std::vector<std::uint8_t> orbit_by_signature;
  std::vector<Q5Orbit> orbits;
  std::uint64_t generator_transitions = 0;
  std::uint64_t map_checksum = 0;
};

struct Profile {
  std::uint64_t low = 0;
  std::uint64_t high = 0;

  auto operator<=>(const Profile&) const = default;
};

struct ProfileWitness {
  Profile profile;
  std::uint64_t witness = 0;
};

struct LocalCensus {
  std::uint64_t tensors_tested = 0;
  std::uint64_t nondegenerate_tensors = 0;
  std::uint64_t minimum_nondegenerate_hyperplanes =
      std::numeric_limits<std::uint64_t>::max();
  std::uint64_t minimum_witness = 0;
  std::array<std::uint64_t, 64> hyperplane_count_histogram{};
  std::array<std::uint64_t, 66> profile_size_histogram{};
  std::vector<ProfileWitness> profiles;
};

[[nodiscard]] std::uint64_t parse_unsigned(const std::string& text) {
  std::size_t parsed = 0;
  const auto value = std::stoull(text, &parsed, 0);
  if (parsed != text.size()) {
    throw std::invalid_argument("invalid unsigned integer: " + text);
  }
  return value;
}

[[nodiscard]] std::vector<std::uint32_t> kernel_basis(
    std::uint32_t q,
    std::uint32_t functional) {
  if (functional == 0 || functional >= (std::uint32_t{1} << q)) {
    throw std::invalid_argument("invalid hyperplane functional");
  }
  const auto pivot = static_cast<std::uint32_t>(std::countr_zero(functional));
  std::vector<std::uint32_t> basis;
  basis.reserve(q - 1);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    if (coordinate == pivot) continue;
    auto vector = std::uint32_t{1} << coordinate;
    if ((functional >> coordinate) & 1U) {
      vector |= std::uint32_t{1} << pivot;
    }
    basis.push_back(vector);
  }
  return basis;
}

[[nodiscard]] PackedAction25 make_q5_action(
    const std::array<std::uint32_t, 5>& basis) {
  std::array<std::uint32_t, 25> images{};
  for (std::uint32_t bit = 0; bit < images.size(); ++bit) {
    images[bit] = static_cast<std::uint32_t>(
        utsp::transform_cubic_tensor_signature(
            5, std::uint64_t{1} << bit, basis));
  }
  PackedAction25 action;
  for (std::uint32_t chunk = 0; chunk < action.chunks.size(); ++chunk) {
    for (std::uint32_t value = 1; value < action.chunks[chunk].size(); ++value) {
      const auto bit = static_cast<std::uint32_t>(std::countr_zero(value));
      action.chunks[chunk][value] =
          action.chunks[chunk][value ^ (std::uint32_t{1} << bit)] ^
          images[5 * chunk + bit];
    }
  }
  return action;
}

[[nodiscard]] Q5Census enumerate_q5_orbits() {
  constexpr std::uint32_t tensor_count = std::uint32_t{1} << 25;
  const std::array<PackedAction25, 2> generators{
      make_q5_action({2, 4, 8, 16, 1}),
      make_q5_action({3, 2, 4, 8, 16}),
  };

  Q5Census result;
  result.orbit_by_signature.assign(tensor_count, kUnassigned);
  std::vector<std::uint32_t> queue;
  queue.reserve(1 << 20);
  for (std::uint32_t start = 0; start < tensor_count; ++start) {
    if (result.orbit_by_signature[start] != kUnassigned) continue;
    if (result.orbits.size() >= kUnassigned) {
      throw std::logic_error("q5 orbit index exceeds packed map");
    }
    const auto orbit_index = static_cast<std::uint8_t>(result.orbits.size());
    queue.clear();
    queue.push_back(start);
    result.orbit_by_signature[start] = orbit_index;
    for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) {
      const auto current = queue[cursor];
      for (const auto& generator : generators) {
        ++result.generator_transitions;
        const auto next = generator(current);
        if (result.orbit_by_signature[next] == kUnassigned) {
          result.orbit_by_signature[next] = orbit_index;
          queue.push_back(next);
        } else if (result.orbit_by_signature[next] != orbit_index) {
          throw std::logic_error("q5 generator connected distinct orbits");
        }
      }
    }
    const auto canonical = utsp::canonicalize_cubic_tensor_direct(5, start);
    const auto radical_dimension = static_cast<std::uint32_t>(
        utsp::cubic_tensor_radical_basis_from_signature(5, start).size());
    result.orbits.push_back(
        {static_cast<std::uint32_t>(canonical.canonical_key),
         start,
         static_cast<std::uint64_t>(queue.size()),
         radical_dimension,
         kUnassigned});
  }

  std::vector<std::size_t> nondegenerate;
  for (std::size_t index = 0; index < result.orbits.size(); ++index) {
    if (result.orbits[index].radical_dimension == 0) {
      nondegenerate.push_back(index);
    }
  }
  std::sort(
      nondegenerate.begin(), nondegenerate.end(),
      [&](std::size_t left, std::size_t right) {
        return result.orbits[left].canonical_key <
               result.orbits[right].canonical_key;
      });
  if (nondegenerate.size() > 128) {
    throw std::logic_error("q5 output profile exceeds 128 packed bits");
  }
  for (std::size_t bit = 0; bit < nondegenerate.size(); ++bit) {
    result.orbits[nondegenerate[bit]].nondegenerate_bit =
        static_cast<std::uint8_t>(bit);
  }

  std::uint64_t checksum = 1469598103934665603ULL;
  for (const auto value : result.orbit_by_signature) {
    checksum ^= value;
    checksum *= 1099511628211ULL;
  }
  result.map_checksum = checksum;
  return result;
}

[[nodiscard]] RestrictionMap make_restriction_map(
    const std::vector<std::uint32_t>& basis) {
  std::array<std::uint32_t, 41> images{};
  for (std::uint32_t bit = 0; bit < images.size(); ++bit) {
    images[bit] = static_cast<std::uint32_t>(
        utsp::restrict_cubic_tensor_signature(
            6, std::uint64_t{1} << bit, basis));
  }
  RestrictionMap map;
  constexpr std::array<std::uint32_t, 3> starts{0, 14, 28};
  constexpr std::array<std::uint32_t, 3> widths{14, 14, 13};
  for (std::size_t chunk = 0; chunk < map.chunks.size(); ++chunk) {
    auto& table = map.chunks[chunk];
    table.resize(std::size_t{1} << widths[chunk]);
    for (std::uint32_t value = 1; value < table.size(); ++value) {
      const auto bit = static_cast<std::uint32_t>(std::countr_zero(value));
      table[value] = table[value ^ (std::uint32_t{1} << bit)] ^
                     images[starts[chunk] + bit];
    }
  }
  return map;
}

[[nodiscard]] std::vector<RestrictionMap> make_q6_hyperplane_maps() {
  std::vector<RestrictionMap> maps;
  maps.reserve(63);
  for (std::uint32_t functional = 1; functional < 64; ++functional) {
    maps.push_back(make_restriction_map(kernel_basis(6, functional)));
  }
  return maps;
}

void set_profile_bit(Profile& profile, std::uint8_t bit) {
  if (bit == kUnassigned) return;
  if (bit < 64) {
    profile.low |= std::uint64_t{1} << bit;
  } else {
    profile.high |= std::uint64_t{1} << (bit - 64);
  }
}

[[nodiscard]] std::uint32_t profile_size(const Profile& profile) {
  return static_cast<std::uint32_t>(
      std::popcount(profile.low) + std::popcount(profile.high));
}

[[nodiscard]] bool subset(const Profile& left, const Profile& right) {
  return (left.low & ~right.low) == 0 && (left.high & ~right.high) == 0;
}

[[nodiscard]] std::pair<Profile, std::uint32_t> tensor_profile(
    std::uint64_t signature,
    const Q5Census& q5,
    const std::vector<RestrictionMap>& restrictions) {
  Profile profile;
  std::uint32_t nondegenerate_hyperplanes = 0;
  for (const auto& restriction : restrictions) {
    const auto restricted = restriction(signature);
    const auto orbit_index = q5.orbit_by_signature[restricted];
    if (orbit_index == kUnassigned || orbit_index >= q5.orbits.size()) {
      throw std::logic_error("q5 restriction has no orbit assignment");
    }
    const auto bit = q5.orbits[orbit_index].nondegenerate_bit;
    if (bit != kUnassigned) {
      ++nondegenerate_hyperplanes;
      set_profile_bit(profile, bit);
    }
  }
  return {profile, nondegenerate_hyperplanes};
}

[[nodiscard]] std::uint64_t gram_repeated_signature(
    std::uint32_t rank,
    bool alternating) {
  if (rank > 6 || (alternating && (rank & 1U)) ||
      (!alternating && rank == 0)) {
    throw std::invalid_argument("invalid symmetric Gram class");
  }
  std::uint64_t signature = 0;
  if (!alternating) {
    for (std::uint32_t coordinate = 0; coordinate < rank; ++coordinate) {
      signature |= std::uint64_t{1} << coordinate;
    }
    return signature;
  }
  std::uint32_t bit = 6;
  for (std::uint32_t left = 0; left < 6; ++left) {
    for (std::uint32_t right = left + 1; right < 6; ++right, ++bit) {
      if (left < rank && (left & 1U) == 0 && right == left + 1) {
        signature |= std::uint64_t{1} << bit;
      }
    }
  }
  return signature;
}

[[nodiscard]] std::vector<std::uint64_t> gram_representatives() {
  std::vector<std::uint64_t> representatives;
  representatives.push_back(gram_repeated_signature(0, true));
  for (std::uint32_t rank = 1; rank <= 6; ++rank) {
    if ((rank & 1U) == 0) {
      representatives.push_back(gram_repeated_signature(rank, true));
    }
    representatives.push_back(gram_repeated_signature(rank, false));
  }
  return representatives;
}

void retain_minimum(
    std::uint64_t value,
    std::uint64_t witness,
    std::uint64_t& minimum,
    std::uint64_t& minimum_witness) {
  if (value < minimum || (value == minimum && witness < minimum_witness)) {
    minimum = value;
    minimum_witness = witness;
  }
}

void sort_and_deduplicate_profiles(std::vector<ProfileWitness>& profiles) {
  std::sort(
      profiles.begin(), profiles.end(),
      [](const ProfileWitness& left, const ProfileWitness& right) {
        if (left.profile != right.profile) return left.profile < right.profile;
        return left.witness < right.witness;
      });
  auto output = profiles.begin();
  for (auto input = profiles.begin(); input != profiles.end();) {
    *output++ = *input;
    const auto profile = input->profile;
    do {
      ++input;
    } while (input != profiles.end() && input->profile == profile);
  }
  profiles.erase(output, profiles.end());
}

[[nodiscard]] std::vector<ProfileWitness> minimal_profiles(
    std::vector<ProfileWitness> profiles) {
  std::sort(
      profiles.begin(), profiles.end(),
      [](const ProfileWitness& left, const ProfileWitness& right) {
        const auto left_size = profile_size(left.profile);
        const auto right_size = profile_size(right.profile);
        if (left_size != right_size) return left_size < right_size;
        if (left.profile != right.profile) return left.profile < right.profile;
        return left.witness < right.witness;
      });
  std::vector<ProfileWitness> retained;
  for (const auto& candidate : profiles) {
    const auto dominated = std::any_of(
        retained.begin(), retained.end(), [&](const ProfileWitness& existing) {
          return subset(existing.profile, candidate.profile);
        });
    if (!dominated) retained.push_back(candidate);
  }
  return retained;
}

template <typename Histogram>
void write_histogram(
    const Histogram& histogram,
    std::uint32_t first_key,
    std::uint32_t last_key) {
  bool wrote = false;
  std::cout << '{';
  for (std::uint32_t key = first_key; key <= last_key; ++key) {
    if (histogram[key] == 0) continue;
    if (wrote) std::cout << ',';
    wrote = true;
    std::cout << '"' << key << "\":" << histogram[key];
  }
  std::cout << '}';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::uint32_t workers = std::max(1U, std::thread::hardware_concurrency());
    std::uint64_t limit = 0;
    for (int index = 1; index < argc; ++index) {
      const std::string argument = argv[index];
      if (argument == "--workers" && index + 1 < argc) {
        workers = static_cast<std::uint32_t>(parse_unsigned(argv[++index]));
      } else if (argument == "--limit" && index + 1 < argc) {
        limit = parse_unsigned(argv[++index]);
      } else if (argument == "--help") {
        std::cout << "usage: q6-hyperplane-output-profile-census "
                     "[--workers N] [--limit N]\n";
        return 0;
      } else {
        throw std::invalid_argument("unknown or incomplete argument: " + argument);
      }
    }
    if (workers == 0) throw std::invalid_argument("workers must be positive");

    const auto start = Clock::now();
    auto q5 = enumerate_q5_orbits();
    if (q5.orbits.size() != 88) {
      throw std::logic_error("q5 orbit census did not find 88 classes");
    }
    const auto q5_nondegenerate = static_cast<std::uint32_t>(std::count_if(
        q5.orbits.begin(), q5.orbits.end(),
        [](const Q5Orbit& orbit) { return orbit.radical_dimension == 0; }));
    if (q5_nondegenerate != 65) {
      throw std::logic_error("q5 orbit census did not find 65 nondegenerate classes");
    }
    const auto restrictions = make_q6_hyperplane_maps();
    const auto grams = gram_representatives();
    if (grams.size() != 10 || restrictions.size() != 63) {
      throw std::logic_error("q6 finite-slice setup has the wrong size");
    }

    constexpr std::uint64_t alternating_count = std::uint64_t{1} << 20;
    const auto exact_total = alternating_count * grams.size();
    const auto total = limit == 0 ? exact_total : std::min(limit, exact_total);
    constexpr std::uint64_t chunk_size = 32;
    std::atomic<std::uint64_t> next{0};
    const auto thread_count = static_cast<std::uint32_t>(
        std::min<std::uint64_t>(workers, (total + chunk_size - 1) / chunk_size));
    std::vector<LocalCensus> local(thread_count);
    for (auto& entry : local) {
      entry.profiles.reserve(static_cast<std::size_t>(total / thread_count + 1));
    }
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (std::uint32_t worker = 0; worker < thread_count; ++worker) {
      threads.emplace_back([&, worker]() {
        auto& result = local[worker];
        while (true) {
          const auto begin = next.fetch_add(chunk_size);
          if (begin >= total) return;
          const auto end = std::min(total, begin + chunk_size);
          for (auto flat = begin; flat < end; ++flat) {
            ++result.tensors_tested;
            const auto gram_index = static_cast<std::size_t>(
                flat / alternating_count);
            const auto alternating = flat % alternating_count;
            const auto signature = grams[gram_index] | (alternating << 21);
            if (!utsp::cubic_tensor_radical_basis_from_signature(6, signature)
                     .empty()) {
              continue;
            }
            ++result.nondegenerate_tensors;
            const auto [profile, hyperplanes] =
                tensor_profile(signature, q5, restrictions);
            if (hyperplanes >= result.hyperplane_count_histogram.size()) {
              throw std::logic_error("q6 hyperplane count exceeds 63");
            }
            ++result.hyperplane_count_histogram[hyperplanes];
            ++result.profile_size_histogram[profile_size(profile)];
            retain_minimum(
                hyperplanes,
                signature,
                result.minimum_nondegenerate_hyperplanes,
                result.minimum_witness);
            result.profiles.push_back({profile, signature});
          }
        }
      });
    }
    for (auto& thread : threads) thread.join();

    LocalCensus combined;
    for (auto& result : local) {
      combined.tensors_tested += result.tensors_tested;
      combined.nondegenerate_tensors += result.nondegenerate_tensors;
      retain_minimum(
          result.minimum_nondegenerate_hyperplanes,
          result.minimum_witness,
          combined.minimum_nondegenerate_hyperplanes,
          combined.minimum_witness);
      for (std::size_t index = 0;
           index < combined.hyperplane_count_histogram.size(); ++index) {
        combined.hyperplane_count_histogram[index] +=
            result.hyperplane_count_histogram[index];
      }
      for (std::size_t index = 0;
           index < combined.profile_size_histogram.size(); ++index) {
        combined.profile_size_histogram[index] +=
            result.profile_size_histogram[index];
      }
      sort_and_deduplicate_profiles(result.profiles);
      combined.profiles.insert(
          combined.profiles.end(),
          result.profiles.begin(),
          result.profiles.end());
    }
    sort_and_deduplicate_profiles(combined.profiles);
    const auto minimal = minimal_profiles(combined.profiles);

    for (const auto& entry : minimal) {
      if (!utsp::cubic_tensor_radical_basis_from_signature(6, entry.witness)
               .empty()) {
        throw std::logic_error("minimal-profile witness is degenerate");
      }
      const auto [profile, hyperplanes] =
          tensor_profile(entry.witness, q5, restrictions);
      if (profile != entry.profile || hyperplanes < 32) {
        throw std::logic_error("minimal-profile witness did not reproduce");
      }
    }

    const auto elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    std::cout << "{\"schema\":\"utsp-q6-q5-output-profile-census-v1\""
              << ",\"status\":\"" << (total == exact_total ? "pass" : "sample")
              << "\",\"workers\":" << thread_count
              << ",\"q5_tensor_count\":" << q5.orbit_by_signature.size()
              << ",\"q5_orbit_count\":" << q5.orbits.size()
              << ",\"q5_nondegenerate_orbit_count\":" << q5_nondegenerate
              << ",\"q5_generator_transitions\":" << q5.generator_transitions
              << ",\"q5_orbit_map_fnv1a64\":\"0x" << std::hex
              << q5.map_checksum << std::dec << "\""
              << ",\"q6_gram_class_count\":" << grams.size()
              << ",\"q6_exact_tensor_count\":" << exact_total
              << ",\"q6_tensors_tested\":" << combined.tensors_tested
              << ",\"q6_nondegenerate_tensors\":"
              << combined.nondegenerate_tensors
              << ",\"minimum_nondegenerate_q5_hyperplanes\":"
              << combined.minimum_nondegenerate_hyperplanes
              << ",\"minimum_hyperplane_witness\":\"0x" << std::hex
              << combined.minimum_witness << std::dec << "\""
              << ",\"distinct_output_profiles\":" << combined.profiles.size()
              << ",\"minimal_output_profiles\":" << minimal.size()
              << ",\"hyperplane_count_histogram\":";
    write_histogram(combined.hyperplane_count_histogram, 0, 63);
    std::cout << ",\"profile_size_histogram\":";
    write_histogram(combined.profile_size_histogram, 0, 65);
    std::cout << ",\"profiles\":[";
    for (std::size_t index = 0; index < minimal.size(); ++index) {
      if (index) std::cout << ',';
      const auto& entry = minimal[index];
      std::cout << "{\"profile_index\":" << index
                << ",\"profile_size\":" << profile_size(entry.profile)
                << ",\"mask_low\":\"0x" << std::hex << std::setw(16)
                << std::setfill('0') << entry.profile.low
                << "\",\"mask_high\":\"0x" << std::setw(16)
                << entry.profile.high
                << "\",\"witness\":\"0x" << entry.witness << std::dec
                << "\",\"q5_output_keys\":[";
      bool wrote_key = false;
      for (const auto& orbit : q5.orbits) {
        const auto bit = orbit.nondegenerate_bit;
        if (bit == kUnassigned) continue;
        const auto present = bit < 64
                                 ? ((entry.profile.low >> bit) & 1U)
                                 : ((entry.profile.high >> (bit - 64)) & 1U);
        if (!present) continue;
        if (wrote_key) std::cout << ',';
        wrote_key = true;
        std::cout << "\"0x" << std::hex << std::setw(16)
                  << std::setfill('0') << orbit.canonical_key << std::dec
                  << "\"";
      }
      std::cout << "]}";
    }
    std::cout << "],\"elapsed_seconds\":" << std::setprecision(12)
              << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "q6 hyperplane output-profile census failed: " << error.what()
              << '\n';
    return 1;
  }
}
