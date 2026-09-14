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

constexpr int kQuotientDimension = 9;
constexpr int kAmbientSize = 1 << kQuotientDimension;
constexpr int kAffineDimension = kQuotientDimension + 1;
constexpr int kQuadraticFeatureCount =
    1 + kQuotientDimension + kQuotientDimension * (kQuotientDimension - 1) / 2;

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
  fs::path output;
  fs::path summary;
  std::uint64_t first = 0;
  std::uint64_t count = 0;
  int multiplicity = -1;
  int expected_weight = 0;
  int target_length = 0;
};

struct BinaryBasis {
  std::array<std::uint64_t, 64> rows{};
  int rank = 0;

  std::uint64_t reduce(std::uint64_t value) const {
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (((value >> pivot) & 1ULL) != 0 && rows[pivot] != 0) {
        value ^= rows[pivot];
      }
    }
    return value;
  }

  bool add(std::uint64_t value) {
    value = reduce(value);
    if (value == 0) return false;
    const int pivot = 63 - __builtin_clzll(value);
    rows[pivot] = value;
    ++rank;
    return true;
  }
};

struct Totals {
  std::array<std::uint64_t, kQuadraticFeatureCount + 1> rank_counts{};
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
    else if (option == "--output") result.output = value;
    else if (option == "--summary") result.summary = value;
    else if (option == "--first") result.first = std::stoull(value);
    else if (option == "--count") result.count = std::stoull(value);
    else if (option == "--multiplicity") result.multiplicity = std::stoi(value);
    else if (option == "--expected-weight") result.expected_weight = std::stoi(value);
    else if (option == "--target-length") result.target_length = std::stoi(value);
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.input.empty() || result.output.empty() || result.summary.empty() ||
      result.count == 0 || (result.multiplicity != 0 && result.multiplicity != 1) ||
      result.expected_weight <= 0 || result.target_length <= 0 ||
      result.expected_weight + 2 * result.multiplicity != result.target_length) {
    throw std::runtime_error(
        "Usage: --input PATH --output PATH --summary PATH --first N --count N "
        "--multiplicity 0|1 --expected-weight N --target-length N");
  }
  return result;
}

std::uint64_t degree_two_vector(int point) {
  std::uint64_t result = 1;
  int feature = 1;
  for (int coordinate = 0; coordinate < kQuotientDimension; ++coordinate) {
    if (((point >> coordinate) & 1) != 0) result |= 1ULL << feature;
    ++feature;
  }
  for (int left = 0; left < kQuotientDimension; ++left) {
    for (int right = left + 1; right < kQuotientDimension; ++right) {
      if (((point >> left) & 1) != 0 && ((point >> right) & 1) != 0) {
        result |= 1ULL << feature;
      }
      ++feature;
    }
  }
  if (feature != kQuadraticFeatureCount) {
    throw std::runtime_error("Quadratic feature count mismatch");
  }
  return result;
}

bool point_is_set(const std::array<std::uint64_t, 8>& support, int point) {
  return ((support[point / 64] >> (point % 64)) & 1ULL) != 0;
}

ProfileRecord profile(const WeightedRecord& source, int expected_weight,
                      int multiplicity,
                      const std::array<std::uint64_t, kAmbientSize>& features) {
  if (source.members == 0) throw std::runtime_error("Source has zero members");
  BinaryBasis quadratic;
  BinaryBasis affine;
  std::array<std::uint64_t, kAffineDimension> affine_moments{};
  int weight = 0;
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!point_is_set(source.support, point)) continue;
    ++weight;
    const auto value = features[point];
    quadratic.add(value);
    affine.add(static_cast<std::uint64_t>((point << 1) | 1));
    affine_moments[0] ^= value;
    for (int coordinate = 0; coordinate < kQuotientDimension; ++coordinate) {
      if (((point >> coordinate) & 1) != 0) {
        affine_moments[coordinate + 1] ^= value;
      }
    }
  }
  if (weight != expected_weight || affine.rank != kAffineDimension) {
    throw std::runtime_error("Source has incorrect weight or affine rank");
  }
  for (const auto moment : affine_moments) {
    if (moment != 0) throw std::runtime_error("Source fails a quadratic lift moment");
  }
  const int lift_dimension = expected_weight - quadratic.rank - kAffineDimension;
  if (lift_dimension < 0 || lift_dimension >= 64) {
    throw std::runtime_error("Invalid lift quotient dimension");
  }
  std::uint64_t compatible = 1;
  if (multiplicity == 1) {
    compatible = 0;
    for (int point = 0; point < kAmbientSize; ++point) {
      if (!point_is_set(source.support, point) &&
          quadratic.reduce(features[point]) == 0) {
        ++compatible;
      }
    }
  }
  if (compatible > std::numeric_limits<std::uint16_t>::max()) {
    throw std::runtime_error("Compatible fibre count does not fit its record");
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
  if (!fs::exists(left) || !fs::exists(right) ||
      fs::file_size(left) != fs::file_size(right)) return false;
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
      throw std::runtime_error("Refusing to replace a different output: " +
                               destination.string());
    }
    fs::remove(temporary);
  } else {
    fs::rename(temporary, destination);
  }
}

void write_summary(const fs::path& path, const Arguments& arguments,
                   const Totals& totals, double elapsed_seconds) {
  const fs::path temporary = path.string() + ".tmp";
  std::ofstream output(temporary, std::ios::trunc);
  output << "{\n"
         << "  \"status\": \"profiled_m10_source_shard_v1\",\n"
         << "  \"target_length\": " << arguments.target_length << ",\n"
         << "  \"first_source_index\": " << arguments.first << ",\n"
         << "  \"source_count\": " << totals.source_count << ",\n"
         << "  \"multiplicity\": " << arguments.multiplicity << ",\n"
         << "  \"expected_weight\": " << arguments.expected_weight << ",\n"
         << "  \"profile_record_bytes\": " << sizeof(ProfileRecord) << ",\n"
         << "  \"contributing_source_count\": "
         << totals.contributing_source_count << ",\n"
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
         << "  \"elapsed_seconds\": " << elapsed_seconds << "\n"
         << "}\n";
  output.close();
  if (!output) throw std::runtime_error("Could not write the profile summary");
  install_immutable(temporary, path);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto input_bytes = fs::file_size(arguments.input);
    if (input_bytes % sizeof(WeightedRecord) != 0 ||
        arguments.first > input_bytes / sizeof(WeightedRecord) ||
        arguments.count > input_bytes / sizeof(WeightedRecord) - arguments.first) {
      throw std::runtime_error("Requested source range is outside the input ledger");
    }
    std::array<std::uint64_t, kAmbientSize> features{};
    for (int point = 0; point < kAmbientSize; ++point) {
      features[point] = degree_two_vector(point);
    }

    if (!arguments.output.parent_path().empty()) {
      fs::create_directories(arguments.output.parent_path());
    }
    if (!arguments.summary.parent_path().empty()) {
      fs::create_directories(arguments.summary.parent_path());
    }
    const fs::path output_temporary = arguments.output.string() + ".tmp";
    std::ifstream input(arguments.input, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(arguments.first * sizeof(WeightedRecord)));
    std::ofstream output(output_temporary, std::ios::binary | std::ios::trunc);
    Totals totals;
    for (std::uint64_t index = 0; index < arguments.count; ++index) {
      WeightedRecord source;
      input.read(reinterpret_cast<char*>(&source), sizeof(source));
      if (!input) throw std::runtime_error("Source ledger is truncated");
      const auto item = profile(source, arguments.expected_weight,
                                arguments.multiplicity, features);
      output.write(reinterpret_cast<const char*>(&item), sizeof(item));
      if (!output) throw std::runtime_error("Could not write a profile record");
      ++totals.source_count;
      ++totals.rank_counts[item.quadratic_rank];
      totals.compatible_fibre_count += item.compatible_fibre_count;
      totals.zero_compatible_source_count += item.compatible_fibre_count == 0;
      const std::uint64_t lifts = 1ULL << item.lift_dimension;
      const std::uint64_t candidates = arguments.multiplicity == 0
                                           ? lifts - 1
                                           : lifts * item.compatible_fibre_count;
      totals.full_rank_candidate_count += candidates;
      totals.contributing_source_count += candidates != 0;
    }
    output.close();
    if (!output || fs::file_size(output_temporary) !=
                       arguments.count * sizeof(ProfileRecord)) {
      throw std::runtime_error("Profile ledger size mismatch");
    }
    install_immutable(output_temporary, arguments.output);
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
