#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr int kQuotientDimension = 8;
constexpr int kAmbientSize = 1 << kQuotientDimension;

struct MaskRecord {
  std::array<std::uint64_t, 4> mask{};
};

struct RepresentativeRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 4> mask{};
  std::uint64_t member_count = 0;
};

struct SourceRecord {
  std::array<std::uint64_t, 4> mask{};
  std::uint8_t source_kind = 0;
  std::uint8_t multiplicity = 0;
  std::uint8_t quadratic_rank = 0;
  std::uint8_t lift_dimension = 0;
  std::uint32_t reserved = 0;
  std::uint64_t source_index = 0;
  std::uint64_t compatible_fibre_set_count = 0;
  std::uint64_t raw_marked_extension_count = 0;
};

static_assert(sizeof(MaskRecord) == 32);
static_assert(sizeof(RepresentativeRecord) == 72);
static_assert(sizeof(SourceRecord) == 64);

struct Arguments {
  fs::path n0_input;
  fs::path n1_input;
  fs::path n2_input;
  fs::path output;
  std::uint64_t n0_count = 0;
  std::uint64_t n1_count = 0;
  std::uint64_t n2_count = 0;
  int threads = 1;
};

struct XorBasis {
  std::array<std::uint64_t, 64> rows{};
  int rank = 0;

  std::uint64_t reduce(std::uint64_t value) const {
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (rows[pivot] != 0 && ((value >> pivot) & 1ULL) != 0) {
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

bool mask_contains(const std::array<std::uint64_t, 4>& mask, int point) {
  return ((mask[point / 64] >> (point % 64)) & 1ULL) != 0;
}

int mask_weight(const std::array<std::uint64_t, 4>& mask) {
  int result = 0;
  for (const auto word : mask) result += __builtin_popcountll(word);
  return result;
}

std::array<std::uint64_t, kAmbientSize> quadratic_vectors() {
  std::array<std::uint64_t, kAmbientSize> result{};
  for (int point = 0; point < kAmbientSize; ++point) {
    std::array<int, kQuotientDimension> coordinates{};
    std::uint64_t value = 1;
    std::uint64_t bit = 1;
    for (int variable = 0; variable < kQuotientDimension; ++variable) {
      coordinates[variable] =
          (point >> (kQuotientDimension - 1 - variable)) & 1;
      bit <<= 1;
      if (coordinates[variable] != 0) value |= bit;
    }
    for (int left = 0; left < kQuotientDimension; ++left) {
      for (int right = left + 1; right < kQuotientDimension; ++right) {
        bit <<= 1;
        if (coordinates[left] != 0 && coordinates[right] != 0) value |= bit;
      }
    }
    result[point] = value;
  }
  return result;
}

bool build_profile(
    const std::array<std::uint64_t, 4>& mask,
    int multiplicity,
    std::uint64_t source_index,
    int expected_weight,
    const std::array<std::uint64_t, kAmbientSize>& quadratic,
    SourceRecord& output) {
  if (mask_weight(mask) != expected_weight) return false;
  XorBasis quadratic_basis;
  XorBasis affine_basis;
  std::array<std::uint64_t, kQuotientDimension + 1> moment_relations{};
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!mask_contains(mask, point)) continue;
    quadratic_basis.add(quadratic[point]);
    affine_basis.add(static_cast<std::uint64_t>(1U | (point << 1)));
    moment_relations[0] ^= quadratic[point];
    for (int bit = 0; bit < kQuotientDimension; ++bit) {
      if (((point >> bit) & 1) != 0) moment_relations[bit + 1] ^= quadratic[point];
    }
  }
  if (affine_basis.rank != kQuotientDimension + 1 ||
      std::any_of(
          moment_relations.begin(), moment_relations.end(),
          [](std::uint64_t value) { return value != 0; })) {
    return false;
  }
  const int lift_dimension =
      expected_weight - quadratic_basis.rank - (kQuotientDimension + 1);
  if (lift_dimension < 0 || lift_dimension >= 63) return false;

  std::array<std::uint64_t, kAmbientSize> external_labels{};
  int external_count = 0;
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!mask_contains(mask, point)) {
      external_labels[external_count++] = quadratic_basis.reduce(quadratic[point]);
    }
  }
  std::sort(
      external_labels.begin(), external_labels.begin() + external_count);
  std::uint64_t compatible = 0;
  if (multiplicity == 0) {
    compatible = 1;
  } else if (multiplicity == 1) {
    compatible = static_cast<std::uint64_t>(std::count(
        external_labels.begin(), external_labels.begin() + external_count, 0ULL));
  } else if (multiplicity == 2) {
    int start = 0;
    while (start < external_count) {
      int end = start + 1;
      while (end < external_count &&
             external_labels[end] == external_labels[start]) {
        ++end;
      }
      const auto count = static_cast<std::uint64_t>(end - start);
      compatible += count * (count - 1) / 2;
      start = end;
    }
  } else {
    return false;
  }
  output = SourceRecord{
      mask,
      static_cast<std::uint8_t>(multiplicity),
      static_cast<std::uint8_t>(multiplicity),
      static_cast<std::uint8_t>(quadratic_basis.rank),
      static_cast<std::uint8_t>(lift_dimension),
      0,
      source_index,
      compatible,
      compatible << lift_dimension,
  };
  return true;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--n0-input") {
      result.n0_input = value;
    } else if (option == "--n0-count") {
      result.n0_count = std::stoull(value);
    } else if (option == "--n1-input") {
      result.n1_input = value;
    } else if (option == "--n1-count") {
      result.n1_count = std::stoull(value);
    } else if (option == "--n2-input") {
      result.n2_input = value;
    } else if (option == "--n2-count") {
      result.n2_count = std::stoull(value);
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--threads") {
      result.threads = std::stoi(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.n0_input.empty() || result.n1_input.empty() ||
      result.n2_input.empty() || result.output.empty() || result.n0_count == 0 ||
      result.n1_count == 0 || result.n2_count == 0 || result.threads < 1) {
    throw std::runtime_error(
        "Required: --n0-input PATH --n0-count N --n1-input PATH --n1-count N "
        "--n2-input PATH --n2-count N --output PATH --threads N");
  }
  return result;
}

std::vector<std::array<std::uint64_t, 4>> load_masks(
    const fs::path& path,
    std::uint64_t expected_count,
    bool representative_format) {
  const auto record_bytes =
      representative_format ? sizeof(RepresentativeRecord) : sizeof(MaskRecord);
  if (!fs::is_regular_file(path) ||
      fs::file_size(path) != expected_count * record_bytes) {
    throw std::runtime_error("A source-mask input has the wrong byte length");
  }
  std::vector<std::array<std::uint64_t, 4>> result(expected_count);
  std::ifstream input(path, std::ios::binary);
  for (std::uint64_t index = 0; index < expected_count; ++index) {
    if (representative_format) {
      RepresentativeRecord record;
      input.read(reinterpret_cast<char*>(&record), sizeof(record));
      result[index] = record.mask;
    } else {
      MaskRecord record;
      input.read(reinterpret_cast<char*>(&record), sizeof(record));
      result[index] = record.mask;
    }
    if (!input) throw std::runtime_error("A source-mask input is truncated");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite the m=9 source ledger");
    }
    const auto n0 = load_masks(arguments.n0_input, arguments.n0_count, true);
    const auto n1 = load_masks(arguments.n1_input, arguments.n1_count, false);
    const auto n2 = load_masks(arguments.n2_input, arguments.n2_count, false);
    std::vector<std::array<std::uint64_t, 4>> masks;
    masks.reserve(n0.size() + n1.size() + n2.size());
    masks.insert(masks.end(), n0.begin(), n0.end());
    masks.insert(masks.end(), n1.begin(), n1.end());
    masks.insert(masks.end(), n2.begin(), n2.end());
    const auto read_finished = std::chrono::steady_clock::now();
    const auto quadratic = quadratic_vectors();
    std::vector<SourceRecord> records(masks.size());
    std::atomic<std::uint64_t> invalid_count{0};
#ifdef _OPENMP
    omp_set_num_threads(arguments.threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(masks.size()); ++index) {
      int multiplicity = 0;
      std::uint64_t source_index = static_cast<std::uint64_t>(index);
      int expected_weight = 54;
      if (index >= static_cast<std::int64_t>(n0.size() + n1.size())) {
        multiplicity = 2;
        source_index -= n0.size() + n1.size();
        expected_weight = 50;
      } else if (index >= static_cast<std::int64_t>(n0.size())) {
        multiplicity = 1;
        source_index -= n0.size();
        expected_weight = 52;
      }
      if (!build_profile(
              masks[static_cast<std::size_t>(index)], multiplicity,
              source_index, expected_weight, quadratic,
              records[static_cast<std::size_t>(index)])) {
        ++invalid_count;
      }
    }
    const auto compute_finished = std::chrono::steady_clock::now();
    if (invalid_count != 0) {
      throw std::runtime_error("One or more source masks failed exact validation");
    }
    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    output.write(
        reinterpret_cast<const char*>(records.data()),
        static_cast<std::streamsize>(records.size() * sizeof(SourceRecord)));
    output.close();
    if (!output || fs::file_size(temporary) != records.size() * sizeof(SourceRecord)) {
      fs::remove(temporary);
      throw std::runtime_error("Could not write the source ledger");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();

    std::array<std::uint64_t, 3> compatible{};
    std::array<std::uint64_t, 3> raw{};
    std::array<std::uint64_t, 3> full_rank{};
    std::array<std::uint64_t, 3> zero_compatible{};
    std::map<std::tuple<int, int, int>, std::uint64_t> profiles;
    for (const auto& record : records) {
      const int multiplicity = record.multiplicity;
      compatible[multiplicity] += record.compatible_fibre_set_count;
      raw[multiplicity] += record.raw_marked_extension_count;
      full_rank[multiplicity] += record.raw_marked_extension_count -
                                 static_cast<std::uint64_t>(multiplicity == 0);
      zero_compatible[multiplicity] += record.compatible_fibre_set_count == 0;
      ++profiles[{multiplicity, record.quadratic_rank, record.lift_dimension}];
    }
    const auto seconds = [](auto begin, auto end) {
      return std::chrono::duration<double>(end - begin).count();
    };
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m09_native_source_ledger_v1\",\n"
              << "  \"source_count\": " << records.size() << ",\n"
              << "  \"source_counts_by_multiplicity\": [" << n0.size() << ", "
              << n1.size() << ", " << n2.size() << "],\n"
              << "  \"compatible_counts_by_multiplicity\": [" << compatible[0]
              << ", " << compatible[1] << ", " << compatible[2] << "],\n"
              << "  \"raw_counts_by_multiplicity\": [" << raw[0] << ", "
              << raw[1] << ", " << raw[2] << "],\n"
              << "  \"full_rank_counts_by_multiplicity\": [" << full_rank[0]
              << ", " << full_rank[1] << ", " << full_rank[2] << "],\n"
              << "  \"zero_compatible_sources_by_multiplicity\": ["
              << zero_compatible[0] << ", " << zero_compatible[1] << ", "
              << zero_compatible[2] << "],\n"
              << "  \"profile_count\": " << profiles.size() << ",\n"
              << "  \"threads\": " << arguments.threads << ",\n"
              << "  \"read_seconds\": " << seconds(started, read_finished) << ",\n"
              << "  \"compute_seconds\": "
              << seconds(read_finished, compute_finished) << ",\n"
              << "  \"write_seconds\": " << seconds(compute_finished, finished)
              << "\n}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
