#include "utsp/tensor_unmarking.hpp"

#include "utsp/tensor_canonical.hpp"
#include "utsp/tensor_extensions.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace utsp {
namespace {

constexpr std::uint32_t kQ5 = 5;
constexpr std::uint32_t kQ6 = 6;
constexpr std::uint32_t kQ7 = 7;
constexpr std::uint32_t kQ6ExtensionBits = 16;
constexpr std::uint32_t kQ7ExtensionBits = 22;
constexpr std::uint32_t kQ6ExtensionStates = 1U << kQ6ExtensionBits;
constexpr std::uint32_t kQ7ExtensionStates = 1U << kQ7ExtensionBits;
constexpr std::uint64_t kGl6Order = 20158709760ULL;
constexpr std::uint64_t kGl7Order = 163849992929280ULL;
constexpr std::uint32_t kSourceBits = 23;
constexpr std::uint32_t kStateBits = 22;
constexpr std::uint64_t kSourceMask = (std::uint64_t{1} << kSourceBits) - 1;
constexpr std::uint64_t kStateMask = (std::uint64_t{1} << kStateBits) - 1;

template <typename Integer>
[[nodiscard]] Integer read_little(std::istream& input) {
  static_assert(std::is_integral_v<Integer>);
  using Unsigned = std::make_unsigned_t<Integer>;
  Unsigned value = 0;
  for (std::size_t byte = 0; byte < sizeof(Integer); ++byte) {
    const int next = input.get();
    if (next == std::char_traits<char>::eof()) {
      throw std::runtime_error("unexpected end of tensor-unmarking input");
    }
    value |= static_cast<Unsigned>(static_cast<unsigned char>(next))
             << (8 * byte);
  }
  return static_cast<Integer>(value);
}

template <typename Integer>
void write_little(std::ostream& output, Integer value) {
  static_assert(std::is_integral_v<Integer>);
  using Unsigned = std::make_unsigned_t<Integer>;
  const auto unsigned_value = static_cast<Unsigned>(value);
  for (std::size_t byte = 0; byte < sizeof(Integer); ++byte) {
    output.put(static_cast<char>((unsigned_value >> (8 * byte)) & 0xffU));
  }
  if (!output) throw std::runtime_error("could not write tensor class map");
}

[[nodiscard]] std::array<char, 8> read_magic(std::istream& input) {
  std::array<char, 8> magic{};
  input.read(magic.data(), magic.size());
  if (!input) throw std::runtime_error("truncated tensor-unmarking magic");
  return magic;
}

[[nodiscard]] std::array<char, 8> magic(const char (&text)[9]) {
  std::array<char, 8> result{};
  std::copy_n(text, 8, result.begin());
  return result;
}

struct Q5ParentRecord {
  std::uint64_t signature = 0;
  std::uint32_t first_q6_marked = 0;
  std::uint32_t q6_marked_count = 0;
};

struct Q6ParentRecord {
  std::uint64_t signature = 0;
  std::uint64_t automorphism_order = 0;
  std::uint64_t first_q7_marked = 0;
  std::uint32_t q7_marked_count = 0;
};

struct TransitionRecord {
  std::uint32_t q5_parent_index = 0;
  std::uint64_t packed_q7_basis = 0;
  std::uint32_t pivot_coordinate = 0;
};

struct Q6MarkedRecord {
  std::uint32_t canonical_q6_parent_index = 0;
  std::uint32_t extension_representative = 0;
  std::uint64_t packed_transport_basis = 0;
};

struct UnmarkAuthority {
  std::uint32_t transition_count = 0;
  std::uint64_t q7_marked_count = 0;
  std::uint32_t expected_q7_unmarked_count = 0;
  std::vector<Q5ParentRecord> q5_parents;
  std::vector<Q6ParentRecord> q6_parents;
  std::vector<TransitionRecord> transitions;
  std::vector<Q6MarkedRecord> q6_marked;
};

[[nodiscard]] UnmarkAuthority read_authority(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open q7 unmark authority");
  if (read_magic(input) != magic("UTSPU7A1")) {
    throw std::runtime_error("q7 unmark authority has wrong magic");
  }
  const auto version = read_little<std::uint32_t>(input);
  const auto q5 = read_little<std::uint32_t>(input);
  const auto q6 = read_little<std::uint32_t>(input);
  const auto q7 = read_little<std::uint32_t>(input);
  const auto q5_count = read_little<std::uint32_t>(input);
  const auto q6_count = read_little<std::uint32_t>(input);
  const auto q6_marked_count = read_little<std::uint32_t>(input);
  const auto transition_count = read_little<std::uint32_t>(input);
  const auto q7_marked_count = read_little<std::uint64_t>(input);
  const auto expected_q7_count = read_little<std::uint64_t>(input);
  if (version != 1 || q5 != kQ5 || q6 != kQ6 || q7 != kQ7 ||
      q5_count != 88 || q6_count != 767 || q6_marked_count != 15332 ||
      transition_count == 0 || transition_count > 126 ||
      q7_marked_count >= (std::uint64_t{1} << kSourceBits) ||
      expected_q7_count > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("q7 unmark authority header is inconsistent");
  }
  UnmarkAuthority authority;
  authority.transition_count = transition_count;
  authority.q7_marked_count = q7_marked_count;
  authority.expected_q7_unmarked_count =
      static_cast<std::uint32_t>(expected_q7_count);
  authority.q5_parents.resize(q5_count);
  for (auto& record : authority.q5_parents) {
    record.signature = read_little<std::uint64_t>(input);
    record.first_q6_marked = read_little<std::uint32_t>(input);
    record.q6_marked_count = read_little<std::uint32_t>(input);
  }
  authority.q6_parents.resize(q6_count);
  for (auto& record : authority.q6_parents) {
    record.signature = read_little<std::uint64_t>(input);
    record.automorphism_order = read_little<std::uint64_t>(input);
    record.first_q7_marked = read_little<std::uint64_t>(input);
    record.q7_marked_count = read_little<std::uint32_t>(input);
    const auto reserved = read_little<std::uint32_t>(input);
    if (reserved != 0) throw std::runtime_error("nonzero q6 authority padding");
  }
  authority.transitions.resize(
      static_cast<std::size_t>(q6_count) * transition_count);
  for (auto& record : authority.transitions) {
    record.q5_parent_index = read_little<std::uint32_t>(input);
    record.packed_q7_basis = read_little<std::uint64_t>(input);
    record.pivot_coordinate = read_little<std::uint32_t>(input);
  }
  authority.q6_marked.resize(q6_marked_count);
  for (auto& record : authority.q6_marked) {
    record.canonical_q6_parent_index = read_little<std::uint32_t>(input);
    record.extension_representative = read_little<std::uint32_t>(input);
    record.packed_transport_basis = read_little<std::uint64_t>(input);
  }
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("q7 unmark authority has trailing bytes");
  }
  std::uint32_t q6_marked_offset = 0;
  for (const auto& record : authority.q5_parents) {
    if (record.first_q6_marked != q6_marked_offset ||
        record.signature >= (std::uint64_t{1} << 25)) {
      throw std::runtime_error("q5 authority records are not consecutive");
    }
    q6_marked_offset += record.q6_marked_count;
  }
  if (q6_marked_offset != q6_marked_count) {
    throw std::runtime_error("q5 authority records do not cover q6 marked data");
  }
  std::uint64_t q7_marked_offset = 0;
  for (const auto& record : authority.q6_parents) {
    if (record.first_q7_marked != q7_marked_offset ||
        record.signature >= (std::uint64_t{1} << 41) ||
        record.automorphism_order == 0 ||
        kGl6Order % record.automorphism_order != 0) {
      throw std::runtime_error("q6 authority records are inconsistent");
    }
    q7_marked_offset += record.q7_marked_count;
  }
  if (q7_marked_offset != q7_marked_count) {
    throw std::runtime_error("q6 authority records do not cover q7 marked data");
  }
  for (std::size_t parent = 0; parent < q6_count; ++parent) {
    for (std::uint32_t transition = 0; transition < transition_count;
         ++transition) {
      const auto& record = authority.transitions[
          parent * transition_count + transition];
      if (record.q5_parent_index >= q5_count ||
          record.pivot_coordinate >= kQ6) {
        throw std::runtime_error("coordinate transition is out of range");
      }
    }
  }
  for (const auto& record : authority.q6_marked) {
    if (record.canonical_q6_parent_index >= q6_count ||
        record.extension_representative >= kQ6ExtensionStates) {
      throw std::runtime_error("q6 marked authority record is out of range");
    }
  }
  return authority;
}

struct MarkedNode {
  std::uint32_t extension_representative = 0;
  std::uint32_t local_orbit_size = 0;
};

struct MarkedLedger {
  std::vector<std::pair<std::uint64_t, std::uint32_t>> parents;
  std::vector<MarkedNode> nodes;
};

[[nodiscard]] MarkedLedger read_marked_ledger(
    const std::string& path, const UnmarkAuthority& authority) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open q7 marked ledger");
  if (read_magic(input) != magic("UTSPML1\0")) {
    throw std::runtime_error("q7 marked ledger has wrong magic");
  }
  const auto version = read_little<std::uint32_t>(input);
  const auto parent_q = read_little<std::uint32_t>(input);
  const auto child_q = read_little<std::uint32_t>(input);
  const auto parent_count = read_little<std::uint32_t>(input);
  const auto node_count = read_little<std::uint64_t>(input);
  if (version != 1 || parent_q != kQ6 || child_q != kQ7 ||
      parent_count != authority.q6_parents.size() ||
      node_count != authority.q7_marked_count) {
    throw std::runtime_error("q7 marked ledger header is inconsistent");
  }
  MarkedLedger ledger;
  ledger.parents.reserve(parent_count);
  for (std::uint32_t parent = 0; parent < parent_count; ++parent) {
    const auto first = read_little<std::uint64_t>(input);
    const auto count = read_little<std::uint32_t>(input);
    const auto reserved = read_little<std::uint32_t>(input);
    if (reserved != 0 || first != authority.q6_parents[parent].first_q7_marked ||
        count != authority.q6_parents[parent].q7_marked_count) {
      throw std::runtime_error("q7 marked parent record disagrees with authority");
    }
    ledger.parents.emplace_back(first, count);
  }
  ledger.nodes.resize(static_cast<std::size_t>(node_count));
  for (auto& node : ledger.nodes) {
    node.extension_representative = read_little<std::uint32_t>(input);
    node.local_orbit_size = read_little<std::uint32_t>(input);
    if (node.extension_representative >= kQ7ExtensionStates ||
        node.local_orbit_size == 0) {
      throw std::runtime_error("q7 marked node is invalid");
    }
  }
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("q7 marked ledger has trailing bytes");
  }
  return ledger;
}

[[nodiscard]] std::string map_path(
    const std::string& directory,
    std::uint32_t parent,
    const std::string& suffix) {
  std::string name = "parent_";
  name.push_back(static_cast<char>('0' + (parent / 100) % 10));
  name.push_back(static_cast<char>('0' + (parent / 10) % 10));
  name.push_back(static_cast<char>('0' + parent % 10));
  name += suffix;
  return (std::filesystem::path(directory) / name).string();
}

struct OwnerMap {
  std::uint32_t parent_q = 0;
  std::uint64_t parent_signature = 0;
  std::uint32_t extension_bits = 0;
  std::uint32_t orbit_count = 0;
  std::vector<std::uint32_t> owner;
};

[[nodiscard]] OwnerMap read_owner_map(
    const std::string& path,
    std::uint32_t expected_q,
    std::uint64_t expected_signature,
    std::uint32_t expected_orbits) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open owner map: " + path);
  if (read_magic(input) != magic("UTSPOM1\0")) {
    throw std::runtime_error("owner map has wrong magic: " + path);
  }
  const auto version = read_little<std::uint32_t>(input);
  OwnerMap result;
  result.parent_q = read_little<std::uint32_t>(input);
  result.parent_signature = read_little<std::uint64_t>(input);
  result.extension_bits = read_little<std::uint32_t>(input);
  const auto state_count = read_little<std::uint32_t>(input);
  result.orbit_count = read_little<std::uint32_t>(input);
  if (version != 1 || result.parent_q != expected_q ||
      result.parent_signature != expected_signature ||
      result.extension_bits != 1 + expected_q + expected_q * (expected_q - 1) / 2 ||
      state_count != (std::uint32_t{1} << result.extension_bits) ||
      result.orbit_count != expected_orbits) {
    throw std::runtime_error("owner map header mismatch: " + path);
  }
  result.owner.resize(state_count);
  for (auto& value : result.owner) {
    value = read_little<std::uint32_t>(input);
    if (value >= result.orbit_count) {
      throw std::runtime_error("owner map value is out of range: " + path);
    }
  }
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("owner map has trailing bytes: " + path);
  }
  return result;
}

struct TransporterMap {
  std::vector<std::uint64_t> basis;
};

[[nodiscard]] TransporterMap read_transporter_map(
    const std::string& path,
    std::uint32_t expected_q,
    std::uint64_t expected_signature,
    std::uint32_t expected_orbits) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open transporter map: " + path);
  if (read_magic(input) != magic("UTSPTM1\0")) {
    throw std::runtime_error("transporter map has wrong magic: " + path);
  }
  const auto version = read_little<std::uint32_t>(input);
  const auto parent_q = read_little<std::uint32_t>(input);
  const auto parent_signature = read_little<std::uint64_t>(input);
  const auto extension_bits = read_little<std::uint32_t>(input);
  const auto state_count = read_little<std::uint32_t>(input);
  const auto orbit_count = read_little<std::uint32_t>(input);
  const auto replay_checks = read_little<std::uint64_t>(input);
  const auto replay_mismatches = read_little<std::uint64_t>(input);
  if (version != 1 || parent_q != expected_q ||
      parent_signature != expected_signature || extension_bits != 16 ||
      state_count != kQ6ExtensionStates || orbit_count != expected_orbits ||
      replay_checks == 0 || replay_mismatches != 0) {
    throw std::runtime_error("transporter map header mismatch: " + path);
  }
  TransporterMap result;
  result.basis.resize(state_count);
  for (auto& basis : result.basis) basis = read_little<std::uint64_t>(input);
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("transporter map has trailing bytes: " + path);
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> unpack_basis(
    std::uint64_t packed, std::uint32_t dimension) {
  const auto mask = (std::uint64_t{1} << dimension) - 1;
  std::vector<std::uint32_t> result(dimension);
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    result[coordinate] = static_cast<std::uint32_t>(
        (packed >> (coordinate * dimension)) & mask);
  }
  return result;
}

[[nodiscard]] std::uint64_t compose_packed_bases(
    std::uint64_t left_packed,
    std::uint64_t right_packed,
    std::uint32_t dimension) {
  const auto mask = (std::uint64_t{1} << dimension) - 1;
  std::uint64_t result = 0;
  for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
    auto coefficients = static_cast<std::uint32_t>(
        (right_packed >> (coordinate * dimension)) & mask);
    std::uint32_t image = 0;
    while (coefficients != 0) {
      const auto bit = static_cast<std::uint32_t>(
          std::countr_zero(coefficients));
      image ^= static_cast<std::uint32_t>(
          (left_packed >> (bit * dimension)) & mask);
      coefficients &= coefficients - 1;
    }
    result |= std::uint64_t{image} << (coordinate * dimension);
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> extend_basis_to_q7(
    std::uint64_t packed_q6_basis) {
  auto basis = unpack_basis(packed_q6_basis, kQ6);
  basis.push_back(1U << kQ6);
  return basis;
}

struct Route {
  std::uint32_t q6_parent = 0;
  std::uint64_t packed_basis = 0;
};

struct PackedAffineMap {
  std::uint32_t offset = 0;
  std::array<std::array<std::uint32_t, 256>, 3> tables{};

  [[nodiscard]] std::uint32_t apply(std::uint32_t value) const {
    return offset ^ tables[0][value & 0xffU] ^
           tables[1][(value >> 8) & 0xffU] ^
           tables[2][(value >> 16) & 0x3fU];
  }
};

[[nodiscard]] PackedAffineMap pack_affine_map(
    std::uint32_t offset,
    const std::array<std::uint32_t, kQ7ExtensionBits>& columns) {
  PackedAffineMap result;
  result.offset = offset;
  for (std::uint32_t chunk = 0; chunk < 3; ++chunk) {
    for (std::uint32_t mask = 1; mask < 256; ++mask) {
      const auto least = mask & (~mask + 1U);
      const auto local = static_cast<std::uint32_t>(std::countr_zero(least));
      const auto bit = 8 * chunk + local;
      result.tables[chunk][mask] = result.tables[chunk][mask ^ least];
      if (bit < columns.size()) result.tables[chunk][mask] ^= columns[bit];
    }
  }
  return result;
}

struct PreparedTransition {
  std::uint32_t q5_parent = 0;
  std::uint64_t packed_q7_basis = 0;
  PackedAffineMap q6_extension;
  PackedAffineMap q7_outer_extension;
};

[[nodiscard]] std::uint32_t transform_outer_extension(
    std::uint32_t state, std::uint64_t packed_q6_basis) {
  std::array<std::uint32_t, kQ6> rows{};
  for (std::uint32_t coordinate = 0; coordinate < kQ6; ++coordinate) {
    if ((state >> (1 + coordinate)) & 1U) {
      rows[coordinate] |= 1U << coordinate;
    }
  }
  std::uint32_t bit = 1 + kQ6;
  for (std::uint32_t first = 0; first < kQ6; ++first) {
    for (std::uint32_t second = first + 1; second < kQ6; ++second, ++bit) {
      if ((state >> bit) & 1U) {
        rows[first] |= 1U << second;
        rows[second] |= 1U << first;
      }
    }
  }
  std::array<std::uint32_t, kQ6> basis{};
  for (std::uint32_t coordinate = 0; coordinate < kQ6; ++coordinate) {
    basis[coordinate] = static_cast<std::uint32_t>(
        (packed_q6_basis >> (coordinate * kQ6)) & 63U);
  }
  std::array<std::uint32_t, kQ6> contracted{};
  for (std::uint32_t coordinate = 0; coordinate < kQ6; ++coordinate) {
    auto vector = basis[coordinate];
    while (vector != 0) {
      const auto row = static_cast<std::uint32_t>(std::countr_zero(vector));
      contracted[coordinate] ^= rows[row];
      vector &= vector - 1;
    }
  }
  std::uint32_t result = state & 1U;
  bit = 1;
  for (std::uint32_t coordinate = 0; coordinate < kQ6; ++coordinate, ++bit) {
    result |= (std::popcount(contracted[coordinate] & basis[coordinate]) & 1U)
              << bit;
  }
  for (std::uint32_t first = 0; first < kQ6; ++first) {
    for (std::uint32_t second = first + 1; second < kQ6; ++second, ++bit) {
      result |= (std::popcount(contracted[first] & basis[second]) & 1U) << bit;
    }
  }
  if (bit != kQ7ExtensionBits) {
    throw std::logic_error("outer extension transform emitted wrong bit count");
  }
  return result;
}

class DisjointSet {
 public:
  explicit DisjointSet(std::size_t size) : parent_(size), size_(size, 1) {
    std::iota(parent_.begin(), parent_.end(), 0U);
  }

  [[nodiscard]] std::uint32_t find(std::uint32_t item) {
    auto root = item;
    while (parent_[root] != root) root = parent_[root];
    while (parent_[item] != item) {
      const auto following = parent_[item];
      parent_[item] = root;
      item = following;
    }
    return root;
  }

  void unite(std::uint32_t left, std::uint32_t right) {
    left = find(left);
    right = find(right);
    if (left == right) return;
    if (size_[left] < size_[right]) std::swap(left, right);
    parent_[right] = left;
    size_[left] += size_[right];
  }

 private:
  std::vector<std::uint32_t> parent_;
  std::vector<std::uint32_t> size_;
};

void write_class_map(
    const std::string& path,
    std::span<const std::uint32_t> class_by_node,
    std::uint32_t class_count) {
  const auto temporary = path + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("could not create tensor class map");
    output.write("UTSPC7M1", 8);
    write_little<std::uint32_t>(output, 1);
    write_little<std::uint32_t>(output, kQ7);
    write_little<std::uint64_t>(output, class_by_node.size());
    write_little<std::uint32_t>(output, class_count);
    write_little<std::uint32_t>(output, 0);
    for (const auto value : class_by_node) {
      write_little<std::uint32_t>(output, value);
    }
  }
  std::filesystem::rename(temporary, path);
}

}  // namespace

TensorUnmarkResult unmark_q7_tensor_orbits(
    const std::string& authority_path,
    const std::string& marked_ledger_path,
    const std::string& q6_map_directory,
    const std::string& q7_owner_map_directory,
    const std::string& class_map_output_path,
    std::uint32_t workers) {
  if (workers == 0) throw std::invalid_argument("worker count must be positive");
  const auto authority = read_authority(authority_path);
  const auto ledger = read_marked_ledger(marked_ledger_path, authority);
  TensorUnmarkResult result;
  result.q5_parent_count = authority.q5_parents.size();
  result.q6_parent_count = authority.q6_parents.size();
  result.q6_marked_orbit_count = authority.q6_marked.size();
  result.transitions_per_marked_orbit = authority.transition_count;
  result.q7_marked_orbit_count = authority.q7_marked_count;
  result.transition_edge_count =
      authority.q7_marked_count * authority.transition_count;
  result.expected_unmarked_orbit_count = authority.expected_q7_unmarked_count;
  result.checks.authority_and_marked_ledger_headers_match = true;

  std::vector<Route> routes(
      authority.q5_parents.size() * kQ6ExtensionStates);
  std::uint64_t route_replay_checks = 0;
  for (std::uint32_t q5_parent = 0;
       q5_parent < authority.q5_parents.size(); ++q5_parent) {
    const auto& parent = authority.q5_parents[q5_parent];
    const auto owner = read_owner_map(
        map_path(q6_map_directory, q5_parent, ".owner.bin"),
        kQ5, parent.signature, parent.q6_marked_count);
    const auto transporter = read_transporter_map(
        map_path(q6_map_directory, q5_parent, ".transporter.bin"),
        kQ5, parent.signature, parent.q6_marked_count);
    for (std::uint32_t local = 0; local < parent.q6_marked_count; ++local) {
      const auto& marked = authority.q6_marked[parent.first_q6_marked + local];
      const auto witness = combine_tensor_parent_and_extension(
          kQ5, parent.signature, marked.extension_representative);
      const auto canonical = transform_cubic_tensor_signature(
          kQ6, witness, unpack_basis(marked.packed_transport_basis, kQ6));
      if (canonical != authority.q6_parents[marked.canonical_q6_parent_index].signature) {
        throw std::runtime_error("q6 marked canonical transport did not replay");
      }
      ++route_replay_checks;
    }
    for (std::uint32_t state = 0; state < kQ6ExtensionStates; ++state) {
      const auto local = owner.owner[state];
      const auto& marked = authority.q6_marked[parent.first_q6_marked + local];
      routes[static_cast<std::size_t>(q5_parent) * kQ6ExtensionStates + state] =
          Route{
              marked.canonical_q6_parent_index,
              compose_packed_bases(
                  transporter.basis[state], marked.packed_transport_basis, kQ6),
          };
    }
  }
  result.checks.every_owner_and_transporter_map_header_matches = true;
  result.checks.every_q6_marked_transport_replays =
      route_replay_checks == authority.q6_marked.size();

  std::vector<PreparedTransition> prepared(authority.transitions.size());
  for (std::uint32_t q6_parent = 0;
       q6_parent < authority.q6_parents.size(); ++q6_parent) {
    const auto parent_signature = authority.q6_parents[q6_parent].signature;
    for (std::uint32_t transition = 0;
         transition < authority.transition_count; ++transition) {
      const auto index =
          static_cast<std::size_t>(q6_parent) * authority.transition_count +
          transition;
      const auto& source = authority.transitions[index];
      const auto basis = unpack_basis(source.packed_q7_basis, kQ7);
      const auto transformed_base = transform_cubic_tensor_signature(
          kQ7,
          combine_tensor_parent_and_extension(kQ6, parent_signature, 0),
          basis);
      const auto transformed_parent =
          restrict_tensor_to_standard_hyperplane(kQ7, transformed_base);
      if (restrict_tensor_to_standard_hyperplane(kQ6, transformed_parent) !=
          authority.q5_parents[source.q5_parent_index].signature) {
        throw std::runtime_error("coordinate transition has wrong q5 restriction");
      }
      const auto q6_offset = tensor_extension_state(kQ5, transformed_parent);
      const auto q7_offset = tensor_extension_state(kQ6, transformed_base);
      std::array<std::uint32_t, kQ7ExtensionBits> q6_columns{};
      std::array<std::uint32_t, kQ7ExtensionBits> q7_columns{};
      for (std::uint32_t bit = 0; bit < kQ7ExtensionBits; ++bit) {
        const auto transformed = transform_cubic_tensor_signature(
            kQ7,
            combine_tensor_parent_and_extension(
                kQ6, parent_signature, std::uint32_t{1} << bit),
            basis);
        q6_columns[bit] =
            tensor_extension_state(
                kQ5, restrict_tensor_to_standard_hyperplane(kQ7, transformed)) ^
            q6_offset;
        q7_columns[bit] = tensor_extension_state(kQ6, transformed) ^ q7_offset;
      }
      prepared[index] = PreparedTransition{
          source.q5_parent_index,
          source.packed_q7_basis,
          pack_affine_map(q6_offset, q6_columns),
          pack_affine_map(q7_offset, q7_columns),
      };
      const auto& parent_nodes = ledger.parents[q6_parent];
      if (parent_nodes.second != 0) {
        const auto sample = ledger.nodes[parent_nodes.first].extension_representative;
        const auto transformed = transform_cubic_tensor_signature(
            kQ7,
            combine_tensor_parent_and_extension(kQ6, parent_signature, sample),
            basis);
        const auto actual_q6 = tensor_extension_state(
            kQ5, restrict_tensor_to_standard_hyperplane(kQ7, transformed));
        const auto actual_q7 = tensor_extension_state(kQ6, transformed);
        if (prepared[index].q6_extension.apply(sample) != actual_q6 ||
            prepared[index].q7_outer_extension.apply(sample) != actual_q7) {
          throw std::runtime_error("transition affine map did not replay");
        }
        ++result.transition_affine_replay_checks;
      }
    }
  }
  result.checks.every_transition_affine_map_replays =
      result.transition_affine_replay_checks == prepared.size();

  std::vector<std::uint64_t> queries(
      static_cast<std::size_t>(result.transition_edge_count));
  std::atomic<std::uint32_t> next_parent{0};
  std::atomic<std::uint64_t> edge_replay_checks{0};
  std::atomic<bool> edge_replay_failed{false};
  const auto worker = [&]() {
    while (true) {
      const auto q6_parent = next_parent.fetch_add(1);
      if (q6_parent >= authority.q6_parents.size()) break;
      const auto& parent = authority.q6_parents[q6_parent];
      const auto [first, count] = ledger.parents[q6_parent];
      for (std::uint32_t local = 0; local < count; ++local) {
        const auto source_node = static_cast<std::uint32_t>(first + local);
        const auto extension = ledger.nodes[source_node].extension_representative;
        for (std::uint32_t transition = 0;
             transition < authority.transition_count; ++transition) {
          const auto& prepared_transition = prepared[
              static_cast<std::size_t>(q6_parent) * authority.transition_count +
              transition];
          const auto q6_state = prepared_transition.q6_extension.apply(extension);
          const auto outer =
              prepared_transition.q7_outer_extension.apply(extension);
          const auto& route = routes[
              static_cast<std::size_t>(prepared_transition.q5_parent) *
                  kQ6ExtensionStates +
              q6_state];
          const auto target_state =
              transform_outer_extension(outer, route.packed_basis);
          const auto packed =
              (std::uint64_t{route.q6_parent} << (kStateBits + kSourceBits)) |
              (std::uint64_t{target_state} << kSourceBits) | source_node;
          queries[static_cast<std::size_t>(source_node) *
                      authority.transition_count +
                  transition] = packed;

          const auto edge_index =
              static_cast<std::uint64_t>(source_node) * authority.transition_count +
              transition;
          const auto replay_stride =
              std::max<std::uint64_t>(1, result.transition_edge_count / 4096);
          if (edge_index % replay_stride == 0) {
            const auto normalized = transform_cubic_tensor_signature(
                kQ7,
                combine_tensor_parent_and_extension(
                    kQ6, parent.signature, extension),
                unpack_basis(prepared_transition.packed_q7_basis, kQ7));
            const auto canonicalized = transform_cubic_tensor_signature(
                kQ7, normalized, extend_basis_to_q7(route.packed_basis));
            if (restrict_tensor_to_standard_hyperplane(kQ7, canonicalized) !=
                    authority.q6_parents[route.q6_parent].signature ||
                tensor_extension_state(kQ6, canonicalized) != target_state) {
              edge_replay_failed.store(true);
            }
            edge_replay_checks.fetch_add(1);
          }
        }
      }
    }
  };
  const auto thread_count = std::min<std::uint32_t>(
      workers, static_cast<std::uint32_t>(authority.q6_parents.size()));
  std::vector<std::thread> threads;
  threads.reserve(thread_count);
  for (std::uint32_t thread = 0; thread < thread_count; ++thread) {
    threads.emplace_back(worker);
  }
  for (auto& thread : threads) thread.join();
  result.transition_edge_replay_checks = edge_replay_checks.load();
  result.checks.sampled_transition_edges_replay_as_full_tensors =
      result.transition_edge_replay_checks != 0 && !edge_replay_failed.load();

  std::sort(queries.begin(), queries.end());
  DisjointSet components(static_cast<std::size_t>(authority.q7_marked_count));
  std::size_t query = 0;
  for (std::uint32_t q6_parent = 0;
       q6_parent < authority.q6_parents.size(); ++q6_parent) {
    const auto& parent = authority.q6_parents[q6_parent];
    const auto owner = read_owner_map(
        map_path(q7_owner_map_directory, q6_parent, ".owner.bin"),
        kQ6, parent.signature, parent.q7_marked_count);
    while (query < queries.size()) {
      const auto packed = queries[query];
      const auto target_parent = static_cast<std::uint32_t>(
          packed >> (kStateBits + kSourceBits));
      if (target_parent != q6_parent) break;
      const auto target_state = static_cast<std::uint32_t>(
          (packed >> kSourceBits) & kStateMask);
      const auto source_node = static_cast<std::uint32_t>(packed & kSourceMask);
      const auto target_node = static_cast<std::uint32_t>(
          parent.first_q7_marked + owner.owner[target_state]);
      components.unite(source_node, target_node);
      ++query;
    }
  }
  if (query != queries.size()) {
    throw std::runtime_error("edge query references an unknown q6 parent");
  }
  result.checks.every_owner_and_transporter_map_header_matches = true;
  queries.clear();
  queries.shrink_to_fit();

  std::vector<std::uint64_t> orbit_mass(authority.q7_marked_count, 0);
  std::vector<std::uint64_t> representative(
      authority.q7_marked_count, std::numeric_limits<std::uint64_t>::max());
  std::vector<std::uint32_t> marked_count(authority.q7_marked_count, 0);
  std::uint64_t total_mass = 0;
  bool complete_parent_mass = true;
  for (std::uint32_t q6_parent = 0;
       q6_parent < authority.q6_parents.size(); ++q6_parent) {
    const auto& parent = authority.q6_parents[q6_parent];
    const auto parent_orbit_size = kGl6Order / parent.automorphism_order;
    const auto [first, count] = ledger.parents[q6_parent];
    std::uint64_t local_mass = 0;
    std::uint32_t previous_extension = 0;
    bool first_extension = true;
    for (std::uint32_t local = 0; local < count; ++local) {
      const auto node_index = static_cast<std::uint32_t>(first + local);
      const auto& node = ledger.nodes[node_index];
      if (!first_extension && node.extension_representative <= previous_extension) {
        throw std::runtime_error("q7 marked representatives are not ordered");
      }
      first_extension = false;
      previous_extension = node.extension_representative;
      local_mass += node.local_orbit_size;
      const auto weight = parent_orbit_size * node.local_orbit_size;
      total_mass += weight;
      const auto root = components.find(node_index);
      orbit_mass[root] += weight;
      representative[root] = std::min(
          representative[root],
          combine_tensor_parent_and_extension(
              kQ6, parent.signature, node.extension_representative));
      ++marked_count[root];
    }
    if (local_mass != kQ7ExtensionStates) complete_parent_mass = false;
  }
  result.checks.every_marked_parent_has_complete_extension_mass =
      complete_parent_mass;

  std::vector<std::uint32_t> roots;
  for (std::uint32_t node = 0; node < authority.q7_marked_count; ++node) {
    if (marked_count[node] != 0) roots.push_back(node);
  }
  std::sort(
      roots.begin(), roots.end(), [&](std::uint32_t left, std::uint32_t right) {
        return representative[left] < representative[right];
      });
  result.checks.component_count_matches_independent_burnside_count =
      roots.size() == authority.expected_q7_unmarked_count;
  result.q7_tensor_mass = total_mass;
  result.checks.unmarked_orbit_mass_is_two_to_63 =
      total_mass == (std::uint64_t{1} << 63);
  result.classes.reserve(roots.size());
  std::vector<std::uint32_t> class_by_root(
      authority.q7_marked_count, std::numeric_limits<std::uint32_t>::max());
  bool sizes_divide = true;
  for (std::uint32_t class_index = 0; class_index < roots.size(); ++class_index) {
    const auto root = roots[class_index];
    if (orbit_mass[root] == 0 || kGl7Order % orbit_mass[root] != 0) {
      sizes_divide = false;
    }
    class_by_root[root] = class_index;
    result.classes.push_back(TensorUnmarkClass{
        class_index,
        representative[root],
        orbit_mass[root],
        orbit_mass[root] == 0 ? 0 : kGl7Order / orbit_mass[root],
        marked_count[root],
    });
  }
  result.checks.every_orbit_size_divides_gl7_order = sizes_divide;
  result.class_by_marked_node.resize(authority.q7_marked_count);
  for (std::uint32_t node = 0; node < authority.q7_marked_count; ++node) {
    const auto root = components.find(node);
    const auto class_index = class_by_root[root];
    if (class_index == std::numeric_limits<std::uint32_t>::max()) {
      throw std::runtime_error("marked node has no unmarked class");
    }
    result.class_by_marked_node[node] = class_index;
  }
  write_class_map(
      class_map_output_path, result.class_by_marked_node, result.classes.size());
  result.checks.class_map_written = true;

  result.pass =
      result.checks.authority_and_marked_ledger_headers_match &&
      result.checks.every_owner_and_transporter_map_header_matches &&
      result.checks.every_q6_marked_transport_replays &&
      result.checks.every_transition_affine_map_replays &&
      result.checks.sampled_transition_edges_replay_as_full_tensors &&
      result.checks.every_marked_parent_has_complete_extension_mass &&
      result.checks.component_count_matches_independent_burnside_count &&
      result.checks.unmarked_orbit_mass_is_two_to_63 &&
      result.checks.every_orbit_size_divides_gl7_order &&
      result.checks.class_map_written;
  return result;
}

}  // namespace utsp
