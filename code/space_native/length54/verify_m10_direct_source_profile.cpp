#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 9;
constexpr int kAmbientSize = 1 << kDimension;
constexpr int kAffineDimension = kDimension + 1;
constexpr int kFeatureCount = 1 + kDimension + kDimension * (kDimension - 1) / 2;

struct WeightedRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 8> support{};
  std::uint64_t members = 0;
};

struct ProfileRecord {
  std::uint8_t quadratic_rank = 0;
  std::uint8_t lift_dimension = 0;
  std::uint16_t compatible_fibre_count = 0;
};

static_assert(sizeof(WeightedRecord) == 104);
static_assert(sizeof(ProfileRecord) == 4);

struct Arguments {
  fs::path input;
  fs::path profiles;
  fs::path summary;
  std::uint64_t first = 0;
  std::uint64_t count = 0;
  int multiplicity = -1;
  int expected_weight = 0;
  int target_length = 0;
};

struct LowPivotBasis {
  std::array<std::uint64_t, 64> rows{};
  int rank = 0;

  std::uint64_t reduce(std::uint64_t value) const {
    for (int pivot = 0; pivot < 64; ++pivot) {
      if (((value >> pivot) & 1ULL) != 0 && rows[pivot] != 0) value ^= rows[pivot];
    }
    return value;
  }

  bool insert(std::uint64_t value) {
    value = reduce(value);
    if (value == 0) return false;
    rows[__builtin_ctzll(value)] = value;
    ++rank;
    return true;
  }
};

struct Totals {
  std::array<std::uint64_t, kFeatureCount + 1> rank_counts{};
  std::uint64_t source_count = 0;
  std::uint64_t contributing_source_count = 0;
  std::uint64_t zero_compatible_source_count = 0;
  unsigned __int128 compatible_fibre_count = 0;
  unsigned __int128 full_rank_candidate_count = 0;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[index];
    if (option == "--input") result.input = value;
    else if (option == "--profiles") result.profiles = value;
    else if (option == "--summary") result.summary = value;
    else if (option == "--first") result.first = std::stoull(value);
    else if (option == "--count") result.count = std::stoull(value);
    else if (option == "--multiplicity") result.multiplicity = std::stoi(value);
    else if (option == "--expected-weight") result.expected_weight = std::stoi(value);
    else if (option == "--target-length") result.target_length = std::stoi(value);
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.input.empty() || result.profiles.empty() || result.summary.empty() ||
      result.count == 0 || (result.multiplicity != 0 && result.multiplicity != 1) ||
      result.expected_weight <= 0 || result.target_length <= 0 ||
      result.expected_weight + 2 * result.multiplicity != result.target_length) {
    throw std::runtime_error(
        "Usage: --input PATH --profiles PATH --summary PATH --first N --count N "
        "--multiplicity 0|1 --expected-weight N --target-length N");
  }
  return result;
}

std::uint64_t independent_feature_vector(int point) {
  std::array<int, kDimension> coordinates{};
  for (int coordinate = 0; coordinate < kDimension; ++coordinate) {
    coordinates[coordinate] = (point >> coordinate) & 1;
  }
  std::uint64_t result = 1;
  int lane = 1;
  for (const int coordinate : coordinates) {
    if (coordinate != 0) result |= 1ULL << lane;
    ++lane;
  }
  for (int right = 1; right < kDimension; ++right) {
    for (int left = 0; left < right; ++left) {
      if (coordinates[left] != 0 && coordinates[right] != 0) {
        result |= 1ULL << lane;
      }
      ++lane;
    }
  }
  if (lane != kFeatureCount) throw std::runtime_error("Feature count mismatch");
  return result;
}

bool contains(const std::array<std::uint64_t, 8>& support, int point) {
  return ((support[point >> 6] >> (point & 63)) & 1ULL) != 0;
}

ProfileRecord recompute(const WeightedRecord& source, int expected_weight,
                        int multiplicity,
                        const std::array<std::uint64_t, kAmbientSize>& features) {
  if (source.members == 0) throw std::runtime_error("Source has zero multiplicity");
  LowPivotBasis quadratic;
  LowPivotBasis affine;
  std::array<std::uint64_t, kAffineDimension> moments{};
  int weight = 0;
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!contains(source.support, point)) continue;
    ++weight;
    const auto feature = features[point];
    quadratic.insert(feature);
    affine.insert(static_cast<std::uint64_t>(1 | (point << 1)));
    moments[0] ^= feature;
    for (int coordinate = 0; coordinate < kDimension; ++coordinate) {
      if (((point >> coordinate) & 1) != 0) moments[coordinate + 1] ^= feature;
    }
  }
  if (weight != expected_weight || affine.rank != kAffineDimension) {
    throw std::runtime_error("Independent source weight or rank check failed");
  }
  if (std::any_of(moments.begin(), moments.end(),
                  [](std::uint64_t value) { return value != 0; })) {
    throw std::runtime_error("Independent source moment check failed");
  }
  const int lift_dimension = expected_weight - quadratic.rank - kAffineDimension;
  if (lift_dimension < 0 || lift_dimension >= 64) {
    throw std::runtime_error("Independent lift dimension is invalid");
  }
  int compatible = 1;
  if (multiplicity == 1) {
    compatible = 0;
    for (int point = 0; point < kAmbientSize; ++point) {
      if (!contains(source.support, point) && quadratic.reduce(features[point]) == 0) {
        ++compatible;
      }
    }
  }
  if (compatible > std::numeric_limits<std::uint16_t>::max()) {
    throw std::runtime_error("Independent compatible count does not fit");
  }
  return ProfileRecord{static_cast<std::uint8_t>(quadratic.rank),
                       static_cast<std::uint8_t>(lift_dimension),
                       static_cast<std::uint16_t>(compatible)};
}

std::string decimal(unsigned __int128 value) {
  if (value == 0) return "0";
  std::string result;
  while (value != 0) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  std::reverse(result.begin(), result.end());
  return result;
}

bool files_equal(const fs::path& left, const fs::path& right) {
  if (!fs::exists(left) || !fs::exists(right) || fs::file_size(left) != fs::file_size(right)) {
    return false;
  }
  std::ifstream a(left, std::ios::binary);
  std::ifstream b(right, std::ios::binary);
  std::array<char, 1 << 20> aa{};
  std::array<char, 1 << 20> bb{};
  while (a && b) {
    a.read(aa.data(), aa.size());
    b.read(bb.data(), bb.size());
    if (a.gcount() != b.gcount() ||
        !std::equal(aa.begin(), aa.begin() + a.gcount(), bb.begin())) return false;
  }
  return true;
}

void install_immutable(const fs::path& temporary, const fs::path& destination) {
  if (fs::exists(destination)) {
    if (!files_equal(temporary, destination)) {
      fs::remove(temporary);
      throw std::runtime_error("Refusing to replace a different verification summary");
    }
    fs::remove(temporary);
  } else {
    fs::rename(temporary, destination);
  }
}

void write_summary(const fs::path& path, const Arguments& arguments,
                   const Totals& totals, double elapsed) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
  const fs::path temporary = path.string() + ".tmp";
  std::ofstream output(temporary, std::ios::trunc);
  output << "{\n"
         << "  \"status\": \"independently_verified_m10_source_profile_shard_v1\",\n"
         << "  \"target_length\": " << arguments.target_length << ",\n"
         << "  \"first_source_index\": " << arguments.first << ",\n"
         << "  \"source_count\": " << totals.source_count << ",\n"
         << "  \"multiplicity\": " << arguments.multiplicity << ",\n"
         << "  \"expected_weight\": " << arguments.expected_weight << ",\n"
         << "  \"profile_record_bytes\": " << sizeof(ProfileRecord) << ",\n"
         << "  \"contributing_source_count\": " << totals.contributing_source_count
         << ",\n"
         << "  \"zero_compatible_source_count\": "
         << totals.zero_compatible_source_count << ",\n"
         << "  \"compatible_fibre_count\": \""
         << decimal(totals.compatible_fibre_count) << "\",\n"
         << "  \"full_rank_marked_extension_count\": \""
         << decimal(totals.full_rank_candidate_count) << "\",\n"
         << "  \"quadratic_rank_counts\": [";
  for (std::size_t rank = 0; rank < totals.rank_counts.size(); ++rank) {
    if (rank != 0) output << ", ";
    output << totals.rank_counts[rank];
  }
  output << "],\n"
         << "  \"elapsed_seconds\": " << elapsed << "\n"
         << "}\n";
  output.close();
  if (!output) throw std::runtime_error("Could not write verification summary");
  install_immutable(temporary, path);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto source_records = fs::file_size(arguments.input) / sizeof(WeightedRecord);
    if (fs::file_size(arguments.input) % sizeof(WeightedRecord) != 0 ||
        arguments.first > source_records || arguments.count > source_records - arguments.first ||
        fs::file_size(arguments.profiles) != arguments.count * sizeof(ProfileRecord)) {
      throw std::runtime_error("Independent verifier input size or range is invalid");
    }
    std::array<std::uint64_t, kAmbientSize> features{};
    for (int point = 0; point < kAmbientSize; ++point) {
      features[point] = independent_feature_vector(point);
    }
    std::ifstream sources(arguments.input, std::ios::binary);
    sources.seekg(static_cast<std::streamoff>(arguments.first * sizeof(WeightedRecord)));
    std::ifstream profiles(arguments.profiles, std::ios::binary);
    Totals totals;
    for (std::uint64_t local_index = 0; local_index < arguments.count; ++local_index) {
      WeightedRecord source;
      ProfileRecord recorded;
      sources.read(reinterpret_cast<char*>(&source), sizeof(source));
      profiles.read(reinterpret_cast<char*>(&recorded), sizeof(recorded));
      if (!sources || !profiles) throw std::runtime_error("Verifier input is truncated");
      const auto expected = recompute(source, arguments.expected_weight,
                                      arguments.multiplicity, features);
      if (recorded.quadratic_rank != expected.quadratic_rank ||
          recorded.lift_dimension != expected.lift_dimension ||
          recorded.compatible_fibre_count != expected.compatible_fibre_count) {
        throw std::runtime_error("Profile mismatch at source index " +
                                 std::to_string(arguments.first + local_index));
      }
      ++totals.source_count;
      ++totals.rank_counts[expected.quadratic_rank];
      totals.compatible_fibre_count += expected.compatible_fibre_count;
      totals.zero_compatible_source_count += expected.compatible_fibre_count == 0;
      const std::uint64_t lifts = 1ULL << expected.lift_dimension;
      const std::uint64_t candidates = arguments.multiplicity == 0
                                           ? lifts - 1
                                           : lifts * expected.compatible_fibre_count;
      totals.full_rank_candidate_count += candidates;
      totals.contributing_source_count += candidates != 0;
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    write_summary(arguments.summary, arguments, totals, elapsed);
    std::cout << "{\"source_count\":" << totals.source_count
              << ",\"full_rank_marked_extension_count\":\""
              << decimal(totals.full_rank_candidate_count)
              << "\",\"elapsed_seconds\":" << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
