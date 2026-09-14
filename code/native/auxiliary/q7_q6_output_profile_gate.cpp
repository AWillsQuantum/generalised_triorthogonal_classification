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
constexpr std::uint64_t kExactExtensionCount = std::uint64_t{1} << 22;

struct Profile {
  std::uint64_t low = 0;
  std::uint64_t high = 0;

  auto operator<=>(const Profile&) const = default;
};

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

struct Q5Orbit {
  std::uint64_t canonical_key = 0;
  std::uint32_t radical_dimension = 0;
  std::uint8_t nondegenerate_bit = kUnassigned;
};

struct Q5Census {
  std::vector<std::uint8_t> orbit_by_signature;
  std::vector<Q5Orbit> orbits;
  std::uint64_t generator_transitions = 0;
  std::uint64_t map_checksum = 0;
};

struct PackedRestrictionMap {
  std::uint32_t source_bits = 0;
  std::uint32_t chunk_bits = 0;
  std::vector<std::vector<std::uint64_t>> chunks;

  [[nodiscard]] std::uint64_t operator()(std::uint64_t signature) const {
    std::uint64_t result = 0;
    std::uint32_t shift = 0;
    for (const auto& chunk : chunks) {
      const auto width = std::min(chunk_bits, source_bits - shift);
      const auto mask = (std::uint64_t{1} << width) - 1;
      result ^= chunk[(signature >> shift) & mask];
      shift += width;
    }
    return result;
  }
};

struct HyperplaneRestriction {
  std::uint32_t functional = 0;
  PackedRestrictionMap map;
};

struct SupportProfile {
  std::vector<std::uint64_t> keys;
};

struct TargetQ6 {
  std::uint64_t canonical_key = 0;
  std::uint64_t standard_signature = 0;
  std::uint64_t embedded_q7_signature = 0;
  Profile q5_profile;
  std::array<std::uint8_t, 88> q5_histogram{};
  std::uint64_t support_profile_mask = 0;
};

struct Counters {
  std::uint64_t candidates = 0;
  std::uint64_t q6_hyperplanes_examined = 0;
  std::uint64_t q6_degenerate_hyperplanes = 0;
  std::uint64_t q6_nondegenerate_hyperplanes = 0;
  std::uint64_t q5_restrictions_examined = 0;
  std::uint64_t q5_profile_rejections = 0;
  std::uint64_t q6_direct_canonicalizations = 0;
  std::uint64_t q6_canonical_key_rejections = 0;
  std::uint64_t support_profile_rejections = 0;
  std::uint64_t q7_radical_tests = 0;
  std::uint64_t q7_degenerate_rejections = 0;

  void add(const Counters& other) {
    candidates += other.candidates;
    q6_hyperplanes_examined += other.q6_hyperplanes_examined;
    q6_degenerate_hyperplanes += other.q6_degenerate_hyperplanes;
    q6_nondegenerate_hyperplanes += other.q6_nondegenerate_hyperplanes;
    q5_restrictions_examined += other.q5_restrictions_examined;
    q5_profile_rejections += other.q5_profile_rejections;
    q6_direct_canonicalizations += other.q6_direct_canonicalizations;
    q6_canonical_key_rejections += other.q6_canonical_key_rejections;
    support_profile_rejections += other.support_profile_rejections;
    q7_radical_tests += other.q7_radical_tests;
    q7_degenerate_rejections += other.q7_degenerate_rejections;
  }
};

struct LocalResult {
  Counters counters;
  std::vector<std::uint64_t> survivors;
  std::vector<std::uint64_t> minimum_witness;
};

[[nodiscard]] std::uint64_t parse_unsigned(const std::string& text) {
  std::size_t parsed = 0;
  const auto value = std::stoull(text, &parsed, 0);
  if (parsed != text.size()) {
    throw std::invalid_argument("invalid unsigned integer: " + text);
  }
  return value;
}

[[nodiscard]] std::vector<std::uint64_t> parse_profile(
    const std::string& text) {
  std::vector<std::uint64_t> keys;
  std::size_t begin = 0;
  while (begin <= text.size()) {
    const auto end = text.find(',', begin);
    const auto token = text.substr(
        begin, end == std::string::npos ? std::string::npos : end - begin);
    if (token.empty()) throw std::invalid_argument("empty profile key");
    keys.push_back(parse_unsigned(token));
    if (end == std::string::npos) break;
    begin = end + 1;
  }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  return keys;
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
        {canonical.canonical_key, radical_dimension, kUnassigned});
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

[[nodiscard]] PackedRestrictionMap make_restriction_map(
    std::uint32_t source_q,
    const std::vector<std::uint32_t>& basis,
    std::uint32_t chunk_bits) {
  const auto source_bits = utsp::cubic_tensor_term_count(source_q);
  std::vector<std::uint64_t> images(source_bits, 0);
  for (std::uint32_t bit = 0; bit < source_bits; ++bit) {
    images[bit] = utsp::restrict_cubic_tensor_signature(
        source_q, std::uint64_t{1} << bit, basis);
  }

  PackedRestrictionMap map;
  map.source_bits = source_bits;
  map.chunk_bits = chunk_bits;
  for (std::uint32_t begin = 0; begin < source_bits; begin += chunk_bits) {
    const auto width = std::min(chunk_bits, source_bits - begin);
    std::vector<std::uint64_t> table(std::size_t{1} << width, 0);
    for (std::uint32_t value = 1; value < table.size(); ++value) {
      const auto bit = static_cast<std::uint32_t>(std::countr_zero(value));
      table[value] = table[value ^ (std::uint32_t{1} << bit)] ^
                     images[begin + bit];
    }
    map.chunks.push_back(std::move(table));
  }
  return map;
}

[[nodiscard]] std::vector<HyperplaneRestriction> make_hyperplane_maps(
    std::uint32_t source_q,
    std::uint32_t chunk_bits) {
  std::vector<HyperplaneRestriction> maps;
  const auto count = (std::uint32_t{1} << source_q) - 1;
  maps.reserve(count);
  for (std::uint32_t functional = 1; functional <= count; ++functional) {
    maps.push_back(
        {functional,
         make_restriction_map(
             source_q, kernel_basis(source_q, functional), chunk_bits)});
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

[[nodiscard]] Profile q5_output_profile(
    std::uint64_t q6_signature,
    const Q5Census& q5,
    const std::vector<HyperplaneRestriction>& q6_to_q5) {
  Profile profile;
  for (const auto& restriction : q6_to_q5) {
    const auto signature = static_cast<std::uint32_t>(restriction.map(q6_signature));
    const auto orbit_index = q5.orbit_by_signature[signature];
    if (orbit_index == kUnassigned || orbit_index >= q5.orbits.size()) {
      throw std::logic_error("q5 restriction has no orbit assignment");
    }
    set_profile_bit(profile, q5.orbits[orbit_index].nondegenerate_bit);
  }
  return profile;
}

[[nodiscard]] std::array<std::uint8_t, 88> q5_output_histogram(
    std::uint64_t q6_signature,
    const Q5Census& q5,
    const std::vector<HyperplaneRestriction>& q6_to_q5) {
  std::array<std::uint8_t, 88> histogram{};
  for (const auto& restriction : q6_to_q5) {
    const auto signature = static_cast<std::uint32_t>(restriction.map(q6_signature));
    const auto orbit_index = q5.orbit_by_signature[signature];
    if (orbit_index == kUnassigned || orbit_index >= q5.orbits.size()) {
      throw std::logic_error("q5 restriction has no orbit assignment");
    }
    ++histogram[orbit_index];
  }
  return histogram;
}

[[nodiscard]] std::uint64_t interleaved_to_standard(
    std::uint32_t q,
    std::uint64_t interleaved) {
  std::vector<std::uint32_t> identity;
  identity.reserve(q);
  for (std::uint32_t coordinate = 0; coordinate < q; ++coordinate) {
    identity.push_back(std::uint32_t{1} << coordinate);
  }
  const auto terms = utsp::cubic_tensor_term_count(q);
  std::uint64_t standard = 0;
  std::uint64_t covered = 0;
  for (std::uint32_t bit = 0; bit < terms; ++bit) {
    const auto image = utsp::interleaved_tensor_key(
        q, std::uint64_t{1} << bit, identity);
    if (std::popcount(image) != 1 || (covered & image) != 0) {
      throw std::logic_error("standard/interleaved conversion is not a permutation");
    }
    covered |= image;
    if ((interleaved & image) != 0) standard |= std::uint64_t{1} << bit;
  }
  const auto limit = (std::uint64_t{1} << terms) - 1;
  if (covered != limit || (interleaved & ~limit) != 0 ||
      utsp::interleaved_tensor_key(q, standard, identity) != interleaved) {
    throw std::logic_error("standard/interleaved conversion did not round trip");
  }
  return standard;
}

[[nodiscard]] std::vector<std::array<std::uint32_t, 3>> tensor_terms(
    std::uint32_t q) {
  std::vector<std::array<std::uint32_t, 3>> terms;
  terms.reserve(utsp::cubic_tensor_term_count(q));
  for (std::uint32_t first = 0; first < q; ++first) {
    terms.push_back({first, first, first});
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      terms.push_back({first, first, second});
    }
  }
  for (std::uint32_t first = 0; first < q; ++first) {
    for (std::uint32_t second = first + 1; second < q; ++second) {
      for (std::uint32_t third = second + 1; third < q; ++third) {
        terms.push_back({first, second, third});
      }
    }
  }
  return terms;
}

[[nodiscard]] std::pair<std::uint64_t, std::vector<std::uint32_t>>
embed_q6_in_q7(std::uint64_t q6_signature) {
  const auto q6_terms = tensor_terms(6);
  const auto q7_terms = tensor_terms(7);
  std::uint64_t embedded = 0;
  std::vector<std::uint32_t> extension_positions;
  for (std::uint32_t q7_bit = 0; q7_bit < q7_terms.size(); ++q7_bit) {
    const auto& term = q7_terms[q7_bit];
    if (term[2] == 6) {
      extension_positions.push_back(q7_bit);
      continue;
    }
    const auto found = std::find(q6_terms.begin(), q6_terms.end(), term);
    if (found == q6_terms.end()) {
      throw std::logic_error("q6 term was not found in q7 embedding");
    }
    const auto q6_bit = static_cast<std::uint32_t>(found - q6_terms.begin());
    if ((q6_signature >> q6_bit) & 1U) {
      embedded |= std::uint64_t{1} << q7_bit;
    }
  }
  if (extension_positions.size() != 22) {
    throw std::logic_error("q6-to-q7 embedding has the wrong codimension");
  }
  return {embedded, extension_positions};
}

[[nodiscard]] std::array<std::array<std::uint64_t, 2048>, 2>
make_extension_tables(const std::vector<std::uint32_t>& positions) {
  if (positions.size() != 22) {
    throw std::invalid_argument("q7 extension table requires 22 positions");
  }
  std::array<std::array<std::uint64_t, 2048>, 2> tables{};
  for (std::uint32_t half = 0; half < 2; ++half) {
    for (std::uint32_t value = 1; value < tables[half].size(); ++value) {
      const auto bit = static_cast<std::uint32_t>(std::countr_zero(value));
      tables[half][value] =
          tables[half][value ^ (std::uint32_t{1} << bit)] |
          (std::uint64_t{1} << positions[11 * half + bit]);
    }
  }
  return tables;
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

[[nodiscard]] int classify_nondegenerate_q6(
    std::uint64_t signature,
    const Q5Census& q5,
    const std::vector<HyperplaneRestriction>& q6_to_q5,
    const std::vector<TargetQ6>& targets,
    bool trust_q6_profile_uniqueness,
    Counters& counters) {
  ++counters.q6_hyperplanes_examined;
  if (!utsp::cubic_tensor_radical_basis_from_signature(6, signature).empty()) {
    ++counters.q6_degenerate_hyperplanes;
    return -2;
  }
  ++counters.q6_nondegenerate_hyperplanes;

  std::uint64_t candidate_mask =
      (std::uint64_t{1} << targets.size()) - 1;
  std::array<std::uint8_t, 88> observed{};
  for (const auto& restriction : q6_to_q5) {
    ++counters.q5_restrictions_examined;
    const auto restricted = static_cast<std::uint32_t>(restriction.map(signature));
    const auto orbit_index = q5.orbit_by_signature[restricted];
    if (orbit_index == kUnassigned || orbit_index >= q5.orbits.size()) {
      throw std::logic_error("q5 restriction has no orbit assignment");
    }
    ++observed[orbit_index];
    for (std::size_t target = 0; target < targets.size(); ++target) {
      if (((candidate_mask >> target) & 1U) != 0 &&
          observed[orbit_index] > targets[target].q5_histogram[orbit_index]) {
        candidate_mask &= ~(std::uint64_t{1} << target);
      }
    }
    if (candidate_mask == 0) {
      ++counters.q5_profile_rejections;
      return -1;
    }
  }
  for (std::size_t target = 0; target < targets.size(); ++target) {
    if (((candidate_mask >> target) & 1U) != 0 &&
        observed != targets[target].q5_histogram) {
      candidate_mask &= ~(std::uint64_t{1} << target);
    }
  }
  if (candidate_mask == 0) {
    ++counters.q5_profile_rejections;
    return -1;
  }

  if (trust_q6_profile_uniqueness) {
    if (std::popcount(candidate_mask) != 1) {
      throw std::logic_error(
          "a trusted q5 profile did not identify one q6 target");
    }
    return static_cast<int>(std::countr_zero(candidate_mask));
  }

  ++counters.q6_direct_canonicalizations;
  const auto key = utsp::canonicalize_cubic_tensor_direct(6, signature).canonical_key;
  for (std::size_t target = 0; target < targets.size(); ++target) {
    if (key == targets[target].canonical_key) {
      return static_cast<int>(target);
    }
  }
  ++counters.q6_canonical_key_rejections;
  return -1;
}

void write_hex(std::uint64_t value) {
  std::cout << "\"0x" << std::hex << std::setw(16) << std::setfill('0')
            << value << std::dec << "\"";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::uint32_t workers = std::max(1U, std::thread::hardware_concurrency());
    std::uint64_t limit = kExactExtensionCount;
    bool q6_profile_census_only = false;
    bool trust_q6_profile_uniqueness = false;
    std::vector<SupportProfile> profiles;
    for (int index = 1; index < argc; ++index) {
      const std::string argument = argv[index];
      if (argument == "--workers" && index + 1 < argc) {
        workers = static_cast<std::uint32_t>(parse_unsigned(argv[++index]));
      } else if (argument == "--limit" && index + 1 < argc) {
        limit = parse_unsigned(argv[++index]);
      } else if (argument == "--profile" && index + 1 < argc) {
        profiles.push_back({parse_profile(argv[++index])});
      } else if (argument == "--q6-profile-census-only") {
        q6_profile_census_only = true;
      } else if (argument == "--trust-certified-q6-profiles") {
        trust_q6_profile_uniqueness = true;
      } else if (argument == "--help") {
        std::cout << "usage: q7-q6-output-profile-gate --profile KEY[,KEY...] "
                     "[--profile ...] [--workers N] [--limit N] "
                     "[--q6-profile-census-only] "
                     "[--trust-certified-q6-profiles]\n";
        return 0;
      } else {
        throw std::invalid_argument("unknown or incomplete argument: " + argument);
      }
    }
    if (workers == 0 || limit == 0 || limit > kExactExtensionCount ||
        profiles.empty() || profiles.size() > 63) {
      throw std::invalid_argument("invalid worker, limit, or profile count");
    }
    for (const auto& profile : profiles) {
      if (profile.keys.empty()) throw std::invalid_argument("empty support profile");
    }

    const auto start = Clock::now();
    auto q5 = enumerate_q5_orbits();
    const auto q5_nondegenerate = static_cast<std::uint32_t>(std::count_if(
        q5.orbits.begin(), q5.orbits.end(),
        [](const Q5Orbit& orbit) { return orbit.radical_dimension == 0; }));
    if (q5.orbits.size() != 88 || q5_nondegenerate != 65 ||
        q5.map_checksum != 0xb30abb236e63ceafULL) {
      throw std::logic_error("q5 orbit census failed its certified invariants");
    }

    const auto q6_to_q5 = make_hyperplane_maps(6, 14);
    if (q6_to_q5.size() != 63) {
      throw std::logic_error("q6 hyperplane map census has the wrong size");
    }

    std::vector<std::uint64_t> unique_keys;
    for (const auto& profile : profiles) {
      unique_keys.insert(unique_keys.end(), profile.keys.begin(), profile.keys.end());
    }
    std::sort(unique_keys.begin(), unique_keys.end());
    unique_keys.erase(std::unique(unique_keys.begin(), unique_keys.end()),
                      unique_keys.end());
    if (unique_keys.size() > 63) {
      throw std::invalid_argument("too many distinct q6 target keys");
    }

    std::vector<TargetQ6> targets;
    std::vector<std::uint32_t> extension_positions;
    for (const auto key : unique_keys) {
      TargetQ6 target;
      target.canonical_key = key;
      target.standard_signature = interleaved_to_standard(6, key);
      const auto canonical = utsp::canonicalize_cubic_tensor_direct(
          6, target.standard_signature);
      if (canonical.canonical_key != key ||
          !utsp::cubic_tensor_radical_basis_from_signature(
               6, target.standard_signature).empty()) {
        throw std::invalid_argument(
            "a supplied q6 key is not a canonical nondegenerate tensor");
      }
      const auto [embedded, positions] =
          embed_q6_in_q7(target.standard_signature);
      target.embedded_q7_signature = embedded;
      if (extension_positions.empty()) {
        extension_positions = positions;
      } else if (extension_positions != positions) {
        throw std::logic_error("q7 extension positions changed by target");
      }
      target.q5_profile =
          q5_output_profile(target.standard_signature, q5, q6_to_q5);
      target.q5_histogram =
          q5_output_histogram(target.standard_signature, q5, q6_to_q5);
      for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
        if (std::binary_search(
                profiles[profile].keys.begin(), profiles[profile].keys.end(), key)) {
          target.support_profile_mask |= std::uint64_t{1} << profile;
        }
      }
      if (target.support_profile_mask == 0) {
        throw std::logic_error("q6 target key belongs to no support profile");
      }
      targets.push_back(target);
    }

    if (trust_q6_profile_uniqueness) {
      for (std::size_t left = 0; left < targets.size(); ++left) {
        for (std::size_t right = left + 1; right < targets.size(); ++right) {
          if (targets[left].q5_profile == targets[right].q5_profile) {
            throw std::invalid_argument(
                "trusted q6 target profiles are not pairwise distinct");
          }
        }
      }
    }

    if (q6_profile_census_only) {
      const auto grams = gram_representatives();
      if (grams.size() != 10) {
        throw std::logic_error("q6 Gram normal-form census has the wrong size");
      }
      constexpr std::uint64_t alternating_count = std::uint64_t{1} << 20;
      const auto total = alternating_count * grams.size();
      constexpr std::uint64_t chunk_size = 32;
      std::atomic<std::uint64_t> next{0};
      const auto thread_count = static_cast<std::uint32_t>(
          std::min<std::uint64_t>(workers, (total + chunk_size - 1) / chunk_size));
      struct CensusLocal {
        std::uint64_t tested = 0;
        std::uint64_t nondegenerate = 0;
        std::vector<std::uint64_t> matches;
        std::vector<std::uint64_t> minimum_witness;
        std::vector<std::uint64_t> canonical_checks;
        std::vector<std::uint64_t> canonical_mismatches;
        std::vector<std::uint64_t> canonical_check_xor;
      };
      std::vector<CensusLocal> local(thread_count);
      for (auto& result : local) {
        result.matches.assign(targets.size(), 0);
        result.minimum_witness.assign(
            targets.size(), std::numeric_limits<std::uint64_t>::max());
        result.canonical_checks.assign(targets.size(), 0);
        result.canonical_mismatches.assign(targets.size(), 0);
        result.canonical_check_xor.assign(targets.size(), 0);
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
              ++result.tested;
              const auto gram = grams[flat / alternating_count];
              const auto alternating = flat % alternating_count;
              const auto signature = gram | (alternating << 21);
              if (!utsp::cubic_tensor_radical_basis_from_signature(6, signature)
                       .empty()) {
                continue;
              }
              ++result.nondegenerate;
              const auto profile =
                  q5_output_profile(signature, q5, q6_to_q5);
              for (std::size_t target = 0; target < targets.size(); ++target) {
                if (profile == targets[target].q5_profile) {
                  ++result.matches[target];
                  result.minimum_witness[target] =
                      std::min(result.minimum_witness[target], signature);
                  const auto canonical =
                      utsp::canonicalize_cubic_tensor_direct(6, signature)
                          .canonical_key;
                  ++result.canonical_checks[target];
                  if (canonical != targets[target].canonical_key) {
                    ++result.canonical_mismatches[target];
                  }
                  result.canonical_check_xor[target] ^=
                      std::rotl(signature * 0x9e3779b97f4a7c15ULL, 17) ^
                      canonical;
                }
              }
            }
          }
        });
      }
      for (auto& thread : threads) thread.join();
      std::uint64_t tested = 0;
      std::uint64_t nondegenerate = 0;
      std::vector<std::uint64_t> matches(targets.size(), 0);
      std::vector<std::uint64_t> minimum_witness(
          targets.size(), std::numeric_limits<std::uint64_t>::max());
      std::vector<std::uint64_t> canonical_checks(targets.size(), 0);
      std::vector<std::uint64_t> canonical_mismatches(targets.size(), 0);
      std::vector<std::uint64_t> canonical_check_xor(targets.size(), 0);
      for (const auto& result : local) {
        tested += result.tested;
        nondegenerate += result.nondegenerate;
        for (std::size_t target = 0; target < targets.size(); ++target) {
          matches[target] += result.matches[target];
          minimum_witness[target] =
              std::min(minimum_witness[target], result.minimum_witness[target]);
          canonical_checks[target] += result.canonical_checks[target];
          canonical_mismatches[target] += result.canonical_mismatches[target];
          canonical_check_xor[target] ^= result.canonical_check_xor[target];
        }
      }
      if (tested != total) {
        throw std::logic_error("q6 Gram normal-form census was incomplete");
      }
      const auto all_profiles_unique = std::all_of(
          targets.begin(), targets.end(), [&](const TargetQ6& left) {
            return std::count_if(
                       targets.begin(), targets.end(),
                       [&](const TargetQ6& right) {
                         return left.q5_profile == right.q5_profile;
                       }) == 1;
          });
      std::uint64_t total_mismatches = 0;
      for (const auto value : canonical_mismatches) {
        total_mismatches += value;
      }
      if (!all_profiles_unique || total_mismatches != 0) {
        throw std::logic_error(
            "q6 target q5 profiles failed the uniqueness certificate");
      }
      const auto elapsed =
          std::chrono::duration<double>(Clock::now() - start).count();
      std::cout << "{\"schema\":\"utsp-q6-target-q5-profile-census-v1\""
                << ",\"status\":\"pass\""
                << ",\"workers\":" << thread_count
                << ",\"q5_tensor_count\":" << q5.orbit_by_signature.size()
                << ",\"q5_orbit_count\":" << q5.orbits.size()
                << ",\"q5_nondegenerate_orbit_count\":" << q5_nondegenerate
                << ",\"q5_orbit_map_fnv1a64\":\"0x" << std::hex
                << q5.map_checksum << std::dec << "\""
                << ",\"q6_gram_class_count\":" << grams.size()
                << ",\"q6_exact_tensor_count\":" << total
                << ",\"q6_tensors_tested\":" << tested
                << ",\"q6_nondegenerate_tensors\":" << nondegenerate
                << ",\"targets\":[";
      for (std::size_t target = 0; target < targets.size(); ++target) {
        if (target) std::cout << ',';
        std::cout << "{\"canonical_key\":";
        write_hex(targets[target].canonical_key);
        std::cout << ",\"q5_profile_low\":";
        write_hex(targets[target].q5_profile.low);
        std::cout << ",\"q5_profile_high\":";
        write_hex(targets[target].q5_profile.high);
        std::cout << ",\"normal_form_matches\":" << matches[target]
                  << ",\"direct_canonical_checks\":"
                  << canonical_checks[target]
                  << ",\"canonical_key_mismatches\":"
                  << canonical_mismatches[target]
                  << ",\"canonical_check_xor\":\"0x" << std::hex
                  << std::setw(16) << std::setfill('0')
                  << canonical_check_xor[target] << std::dec << "\""
                  << ",\"minimum_standard_signature\":";
        if (minimum_witness[target] ==
            std::numeric_limits<std::uint64_t>::max()) {
          std::cout << "null";
        } else {
          write_hex(minimum_witness[target]);
        }
        std::cout << '}';
      }
      std::cout << "],\"all_target_profiles_distinct\":true"
                << ",\"all_matching_normal_forms_have_claimed_key\":true"
                << ",\"elapsed_seconds\":" << std::setprecision(12)
                << elapsed << "}\n";
      return 0;
    }

    const auto q7_to_q6 = make_hyperplane_maps(7, 9);
    if (q7_to_q6.size() != 127) {
      throw std::logic_error("q7 hyperplane map census has the wrong size");
    }

    const auto anchor = std::find_if(
        q7_to_q6.begin(), q7_to_q6.end(),
        [](const HyperplaneRestriction& restriction) {
          return restriction.functional == 64;
        });
    if (anchor == q7_to_q6.end()) {
      throw std::logic_error("coordinate q6 anchor hyperplane is absent");
    }
    for (const auto& target : targets) {
      if (anchor->map(target.embedded_q7_signature) != target.standard_signature) {
        throw std::logic_error("q7 embedding did not preserve its q6 anchor");
      }
    }
    const auto extension_tables = make_extension_tables(extension_positions);

    const auto total = static_cast<std::uint64_t>(targets.size()) * limit;
    constexpr std::uint64_t chunk_size = 64;
    std::atomic<std::uint64_t> next{0};
    const auto thread_count = static_cast<std::uint32_t>(
        std::min<std::uint64_t>(workers, (total + chunk_size - 1) / chunk_size));
    std::vector<LocalResult> local(thread_count);
    for (auto& result : local) {
      result.survivors.assign(profiles.size(), 0);
      result.minimum_witness.assign(
          profiles.size(), std::numeric_limits<std::uint64_t>::max());
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
            ++result.counters.candidates;
            const auto target_index = static_cast<std::size_t>(flat / limit);
            const auto extension = flat % limit;
            const auto signature =
                targets[target_index].embedded_q7_signature ^
                extension_tables[0][extension & 2047U] ^
                extension_tables[1][(extension >> 11) & 2047U];
            auto profile_mask = targets[target_index].support_profile_mask;
            for (const auto& restriction : q7_to_q6) {
              if (restriction.functional == 64) continue;
              const auto q6_signature = restriction.map(signature);
              const auto classification = classify_nondegenerate_q6(
                  q6_signature,
                  q5,
                  q6_to_q5,
                  targets,
                  trust_q6_profile_uniqueness,
                  result.counters);
              if (classification == -2) continue;
              if (classification < 0) {
                profile_mask = 0;
                break;
              }
              profile_mask &=
                  targets[static_cast<std::size_t>(classification)]
                      .support_profile_mask;
              if (profile_mask == 0) {
                ++result.counters.support_profile_rejections;
                break;
              }
            }
            if (profile_mask == 0) continue;
            ++result.counters.q7_radical_tests;
            if (!utsp::cubic_tensor_radical_basis_from_signature(7, signature)
                     .empty()) {
              ++result.counters.q7_degenerate_rejections;
              continue;
            }
            while (profile_mask != 0) {
              const auto profile = static_cast<std::size_t>(
                  std::countr_zero(profile_mask));
              ++result.survivors[profile];
              result.minimum_witness[profile] =
                  std::min(result.minimum_witness[profile], signature);
              profile_mask &= profile_mask - 1;
            }
          }
        }
      });
    }
    for (auto& thread : threads) thread.join();

    Counters counters;
    std::vector<std::uint64_t> survivors(profiles.size(), 0);
    std::vector<std::uint64_t> minimum_witness(
        profiles.size(), std::numeric_limits<std::uint64_t>::max());
    for (const auto& result : local) {
      counters.add(result.counters);
      for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
        survivors[profile] += result.survivors[profile];
        minimum_witness[profile] =
            std::min(minimum_witness[profile], result.minimum_witness[profile]);
      }
    }
    if (counters.candidates != total) {
      throw std::logic_error("q7 extension enumeration was incomplete");
    }

    const auto elapsed = std::chrono::duration<double>(Clock::now() - start).count();
    std::cout << "{\"schema\":\"utsp-q7-q6-output-profile-gate-v1\""
              << ",\"status\":\""
              << (limit == kExactExtensionCount ? "pass" : "sample") << "\""
              << ",\"workers\":" << thread_count
              << ",\"q5_tensor_count\":" << q5.orbit_by_signature.size()
              << ",\"q5_orbit_count\":" << q5.orbits.size()
              << ",\"q5_nondegenerate_orbit_count\":" << q5_nondegenerate
              << ",\"q5_generator_transitions\":" << q5.generator_transitions
              << ",\"q5_orbit_map_fnv1a64\":\"0x" << std::hex
              << q5.map_checksum << std::dec << "\""
              << ",\"q6_target_keys\":[";
    for (std::size_t index = 0; index < targets.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << "{\"canonical_key\":";
      write_hex(targets[index].canonical_key);
      std::cout << ",\"standard_signature\":";
      write_hex(targets[index].standard_signature);
      std::cout << ",\"q5_profile_low\":";
      write_hex(targets[index].q5_profile.low);
      std::cout << ",\"q5_profile_high\":";
      write_hex(targets[index].q5_profile.high);
      std::cout << ",\"q5_histogram\":[";
      for (std::size_t orbit = 0;
           orbit < targets[index].q5_histogram.size(); ++orbit) {
        if (orbit) std::cout << ',';
        std::cout << static_cast<std::uint32_t>(
            targets[index].q5_histogram[orbit]);
      }
      std::cout << ']';
      std::cout << '}';
    }
    std::cout << "],\"support_profiles\":[";
    for (std::size_t profile = 0; profile < profiles.size(); ++profile) {
      if (profile) std::cout << ',';
      std::cout << "{\"profile_index\":" << profile << ",\"q6_keys\":[";
      for (std::size_t key = 0; key < profiles[profile].keys.size(); ++key) {
        if (key) std::cout << ',';
        write_hex(profiles[profile].keys[key]);
      }
      std::cout << "],\"surviving_anchored_extensions\":" << survivors[profile]
                << ",\"minimum_q7_standard_signature\":";
      if (minimum_witness[profile] == std::numeric_limits<std::uint64_t>::max()) {
        std::cout << "null";
      } else {
        write_hex(minimum_witness[profile]);
      }
      std::cout << '}';
    }
    std::cout << "]"
              << ",\"q6_target_key_count\":" << targets.size()
              << ",\"trusted_certified_q6_profile_uniqueness\":"
              << (trust_q6_profile_uniqueness ? "true" : "false")
              << ",\"extensions_per_anchor\":" << limit
              << ",\"exact_extensions_per_anchor\":" << kExactExtensionCount
              << ",\"candidates_tested\":" << counters.candidates
              << ",\"q6_hyperplanes_examined\":"
              << counters.q6_hyperplanes_examined
              << ",\"q6_degenerate_hyperplanes\":"
              << counters.q6_degenerate_hyperplanes
              << ",\"q6_nondegenerate_hyperplanes\":"
              << counters.q6_nondegenerate_hyperplanes
              << ",\"q5_restrictions_examined\":"
              << counters.q5_restrictions_examined
              << ",\"q5_profile_rejections\":"
              << counters.q5_profile_rejections
              << ",\"q6_direct_canonicalizations\":"
              << counters.q6_direct_canonicalizations
              << ",\"q6_canonical_key_rejections\":"
              << counters.q6_canonical_key_rejections
              << ",\"support_profile_rejections\":"
              << counters.support_profile_rejections
              << ",\"q7_radical_tests\":" << counters.q7_radical_tests
              << ",\"q7_degenerate_rejections\":"
              << counters.q7_degenerate_rejections
              << ",\"elapsed_seconds\":" << std::setprecision(12)
              << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "q7/q6 output-profile gate failed: " << error.what() << '\n';
    return 1;
  }
}
