#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 11;
constexpr int kAmbientSize = 1 << kDimension;
constexpr int kAffineDimension = kDimension + 1;
constexpr int kTargetLength = 54;
constexpr int kFeatureCount =
    1 + kDimension + kDimension * (kDimension - 1) / 2;

struct WeightedRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 32> support{};
  std::uint64_t members = 0;
};

struct ProfileRecord {
  std::uint8_t quadratic_rank = 0;
  std::uint8_t lift_dimension = 0;
  std::uint16_t compatible_fibre_count = 1;
};

static_assert(sizeof(WeightedRecord) == 296);
static_assert(sizeof(ProfileRecord) == 4);

using Feature = unsigned __int128;

struct Arguments {
  fs::path input;
  fs::path profiles;
  fs::path summary;
  std::uint64_t first = 0;
  std::uint64_t count = 0;
};

struct LowPivotBasis {
  std::array<Feature, kFeatureCount> rows{};
  int rank = 0;

  Feature reduce(Feature value) const {
    for (int pivot = 0; pivot < kFeatureCount; ++pivot) {
      if (((value >> pivot) & Feature{1}) != 0 && rows[pivot] != 0) {
        value ^= rows[pivot];
      }
    }
    return value;
  }

  bool insert(Feature value) {
    value = reduce(value);
    if (value == 0) return false;
    for (int pivot = 0; pivot < kFeatureCount; ++pivot) {
      if (((value >> pivot) & Feature{1}) != 0) {
        rows[pivot] = value;
        break;
      }
    }
    ++rank;
    return true;
  }
};

struct Totals {
  std::array<std::uint64_t, kFeatureCount + 1> rank_counts{};
  std::uint64_t source_count = 0;
  std::uint64_t contributing_source_count = 0;
  unsigned __int128 raw_marked_extension_count = 0;
  unsigned __int128 full_rank_marked_extension_count = 0;
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
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.input.empty() || result.profiles.empty() ||
      result.summary.empty() || result.count == 0) {
    throw std::runtime_error(
        "Usage: --input PATH --profiles PATH --summary PATH --first N --count N");
  }
  return result;
}

Feature independent_feature_vector(int point) {
  std::array<int, kDimension> coordinates{};
  for (int coordinate = 0; coordinate < kDimension; ++coordinate) {
    coordinates[coordinate] = (point >> coordinate) & 1;
  }
  Feature result = 1;
  int feature = 1;
  for (const int coordinate : coordinates) {
    if (coordinate != 0) result |= Feature{1} << feature;
    ++feature;
  }
  for (int right = 1; right < kDimension; ++right) {
    for (int left = 0; left < right; ++left) {
      if (coordinates[left] != 0 && coordinates[right] != 0) {
        result |= Feature{1} << feature;
      }
      ++feature;
    }
  }
  if (feature != kFeatureCount) {
    throw std::runtime_error("Independent feature count mismatch");
  }
  return result;
}

ProfileRecord recompute(
    const WeightedRecord& source,
    const std::array<Feature, kAmbientSize>& features) {
  if (source.members == 0) throw std::runtime_error("Source has zero members");
  LowPivotBasis quadratic;
  LowPivotBasis affine;
  std::array<Feature, kAffineDimension> affine_moments{};
  int weight = 0;
  for (int word_index = 0; word_index < 32; ++word_index) {
    std::uint64_t remaining = source.support[word_index];
    while (remaining != 0) {
      const int bit = __builtin_ctzll(remaining);
      const int point = 64 * word_index + bit;
      remaining ^= 1ULL << bit;
      ++weight;
      const auto feature = features[point];
      quadratic.insert(feature);
      affine.insert(static_cast<Feature>(1 | (point << 1)));
      affine_moments[0] ^= feature;
      for (int coordinate = 0; coordinate < kDimension; ++coordinate) {
        if (((point >> coordinate) & 1) != 0) {
          affine_moments[coordinate + 1] ^= feature;
        }
      }
    }
  }
  if (weight != kTargetLength || affine.rank != kAffineDimension) {
    throw std::runtime_error("Independent source weight or rank check failed");
  }
  if (std::any_of(affine_moments.begin(), affine_moments.end(),
                  [](Feature value) { return value != 0; })) {
    throw std::runtime_error("Independent source moment check failed");
  }
  const int lift_dimension = kTargetLength - quadratic.rank - kAffineDimension;
  if (lift_dimension < 0 || lift_dimension >= 64) {
    throw std::runtime_error("Independent lift dimension is invalid");
  }
  return ProfileRecord{static_cast<std::uint8_t>(quadratic.rank),
                       static_cast<std::uint8_t>(lift_dimension), 1};
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
  if (!fs::exists(left) || !fs::exists(right) ||
      fs::file_size(left) != fs::file_size(right)) {
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
        !std::equal(aa.begin(), aa.begin() + a.gcount(), bb.begin())) {
      return false;
    }
  }
  return true;
}

void install_immutable(const fs::path& temporary, const fs::path& destination) {
  if (fs::exists(destination)) {
    if (!files_equal(temporary, destination)) {
      fs::remove(temporary);
      throw std::runtime_error("Refusing to replace a different summary");
    }
    fs::remove(temporary);
  } else {
    fs::rename(temporary, destination);
  }
}

void write_summary(const fs::path& path, const Arguments& arguments,
                   const Totals& totals, double elapsed_seconds) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
  const fs::path temporary = path.string() + ".tmp";
  std::ofstream output(temporary, std::ios::trunc);
  output << "{\n"
         << "  \"status\": \"independently_verified_length54_m12_source_profile_shard_v1\",\n"
         << "  \"target_length\": 54,\n"
         << "  \"target_affine_dimension\": 12,\n"
         << "  \"first_source_index\": " << arguments.first << ",\n"
         << "  \"source_count\": " << totals.source_count << ",\n"
         << "  \"profile_record_bytes\": " << sizeof(ProfileRecord) << ",\n"
         << "  \"contributing_source_count\": "
         << totals.contributing_source_count << ",\n"
         << "  \"raw_marked_extension_count\": \""
         << decimal(totals.raw_marked_extension_count) << "\",\n"
         << "  \"full_rank_marked_extension_count\": \""
         << decimal(totals.full_rank_marked_extension_count) << "\",\n"
         << "  \"quadratic_rank_counts\": [";
  for (std::size_t rank = 0; rank < totals.rank_counts.size(); ++rank) {
    if (rank != 0) output << ", ";
    output << totals.rank_counts[rank];
  }
  output << "],\n"
         << "  \"elapsed_seconds\": " << elapsed_seconds << "\n"
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
    const auto input_bytes = fs::file_size(arguments.input);
    const auto source_records = input_bytes / sizeof(WeightedRecord);
    if (input_bytes % sizeof(WeightedRecord) != 0 ||
        arguments.first > source_records ||
        arguments.count > source_records - arguments.first ||
        fs::file_size(arguments.profiles) !=
            arguments.count * sizeof(ProfileRecord)) {
      throw std::runtime_error("Verifier input size or range is invalid");
    }
    std::array<Feature, kAmbientSize> features{};
    for (int point = 0; point < kAmbientSize; ++point) {
      features[point] = independent_feature_vector(point);
    }
    std::ifstream sources(arguments.input, std::ios::binary);
    sources.seekg(static_cast<std::streamoff>(arguments.first *
                                              sizeof(WeightedRecord)));
    std::ifstream profiles(arguments.profiles, std::ios::binary);
    Totals totals;
    for (std::uint64_t local_index = 0; local_index < arguments.count;
         ++local_index) {
      WeightedRecord source;
      ProfileRecord recorded;
      sources.read(reinterpret_cast<char*>(&source), sizeof(source));
      profiles.read(reinterpret_cast<char*>(&recorded), sizeof(recorded));
      if (!sources || !profiles) {
        throw std::runtime_error("Verifier input is truncated");
      }
      const auto expected = recompute(source, features);
      if (recorded.quadratic_rank != expected.quadratic_rank ||
          recorded.lift_dimension != expected.lift_dimension ||
          recorded.compatible_fibre_count != 1) {
        throw std::runtime_error("Profile mismatch at source index " +
                                 std::to_string(arguments.first + local_index));
      }
      ++totals.source_count;
      ++totals.rank_counts[expected.quadratic_rank];
      const std::uint64_t raw = 1ULL << expected.lift_dimension;
      totals.raw_marked_extension_count += raw;
      totals.full_rank_marked_extension_count += raw - 1;
      totals.contributing_source_count += expected.lift_dimension != 0;
    }
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    write_summary(arguments.summary, arguments, totals, elapsed);
    std::cout << "{\"source_count\":" << totals.source_count
              << ",\"full_rank_marked_extension_count\":\""
              << decimal(totals.full_rank_marked_extension_count)
              << "\",\"elapsed_seconds\":" << elapsed << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
