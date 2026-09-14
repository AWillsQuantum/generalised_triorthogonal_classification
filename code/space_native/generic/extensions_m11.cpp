#include "configuration.hpp"
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
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr int kQuotientDimension = 10;
constexpr int kTargetDimension = 11;
constexpr int kAmbientSize = 1 << kTargetDimension;
constexpr int kQuotientSize = 1 << kQuotientDimension;
constexpr int kSupportSize = UTSP_LENGTH;
constexpr int kProfileBins = kSupportSize / 2 + 1;

struct SourceRecord {
  std::array<std::uint64_t, 16> words{};
  std::uint8_t source_kind = 0;
  std::uint8_t multiplicity = 0;
  std::uint8_t quadratic_rank = 0;
  std::uint8_t lift_dimension = 0;
  std::array<std::uint8_t, 4> padding{};
  std::uint64_t source_index = 0;
  std::uint64_t compatible_fibre_count = 0;
  std::uint64_t raw_marked_extension_count = 0;
};

struct Mask {
  std::array<std::uint64_t, 32> words{};
};

struct SignatureRecord {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

static_assert(sizeof(SourceRecord) == 160);
static_assert(sizeof(Mask) == 256);
static_assert(sizeof(SignatureRecord) == 288);

struct Arguments {
  fs::path input;
  fs::path masks_output;
  fs::path output;
  std::uint64_t source_start = 0;
  std::uint64_t source_count = 0;
  int threads = 1;
};

struct SolverRow {
  std::uint64_t value = 0;
  std::uint64_t coefficients = 0;
};

class ColumnSolver {
 public:
  void add(std::uint64_t column, int index) {
    std::uint64_t value = column;
    std::uint64_t coefficients = 1ULL << index;
    reduce_with_coefficients(value, coefficients);
    if (value == 0) {
      relations_.push_back(coefficients);
      return;
    }
    const int pivot = 63 - __builtin_clzll(value);
    rows_[pivot] = SolverRow{value, coefficients};
    ++rank_;
  }

  std::uint64_t reduce(std::uint64_t value) const {
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (((value >> pivot) & 1ULL) != 0 && rows_[pivot].value != 0) {
        value ^= rows_[pivot].value;
      }
    }
    return value;
  }

  bool solve(std::uint64_t target, std::uint64_t& coefficients) const {
    coefficients = 0;
    reduce_with_coefficients(target, coefficients);
    return target == 0;
  }

  std::uint64_t evaluate(std::uint64_t coefficients,
                         const std::vector<std::uint64_t>& columns) const {
    std::uint64_t result = 0;
    while (coefficients != 0) {
      const std::uint64_t bit = coefficients & (~coefficients + 1);
      result ^= columns[__builtin_ctzll(bit)];
      coefficients ^= bit;
    }
    return result;
  }

  int rank() const { return rank_; }
  const std::vector<std::uint64_t>& relations() const { return relations_; }

 private:
  void reduce_with_coefficients(std::uint64_t& value,
                                std::uint64_t& coefficients) const {
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (((value >> pivot) & 1ULL) != 0 && rows_[pivot].value != 0) {
        value ^= rows_[pivot].value;
        coefficients ^= rows_[pivot].coefficients;
      }
    }
  }

  std::array<SolverRow, 64> rows_{};
  std::vector<std::uint64_t> relations_;
  int rank_ = 0;
};

class CoefficientBasis {
 public:
  bool add(std::uint64_t value) {
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (((value >> pivot) & 1ULL) != 0 && rows_[pivot] != 0) {
        value ^= rows_[pivot];
      }
    }
    if (value == 0) return false;
    rows_[63 - __builtin_clzll(value)] = value;
    ++rank_;
    return true;
  }
  int rank() const { return rank_; }

 private:
  std::array<std::uint64_t, 64> rows_{};
  int rank_ = 0;
};

class ProfileHash {
 public:
  ProfileHash()
      : state_{0xcbf29ce484222325ULL, 0x84222325cbf29ce4ULL,
               0x9e3779b97f4a7c15ULL, 0x6a09e667f3bcc909ULL} {}

  void add_byte(std::uint8_t value) {
    constexpr std::array<std::uint64_t, 4> primes{
        0x00000100000001b3ULL, 0x9e3779b185ebca87ULL,
        0xc2b2ae3d27d4eb4fULL, 0x165667b19e3779f9ULL};
    for (std::size_t index = 0; index < state_.size(); ++index) {
      state_[index] ^=
          static_cast<std::uint64_t>(value) + 0x9dU * index;
      state_[index] *= primes[index];
      state_[index] ^= state_[index] >> (29 + index);
    }
  }

  void add_u16(std::uint16_t value) {
    add_byte(static_cast<std::uint8_t>(value & 0xffU));
    add_byte(static_cast<std::uint8_t>(value >> 8));
  }

  std::array<std::uint64_t, 4> finish() const { return state_; }

 private:
  std::array<std::uint64_t, 4> state_;
};

bool point_is_set(const std::array<std::uint64_t, 16>& words, int point) {
  return ((words[point / 64] >> (point % 64)) & 1ULL) != 0;
}

bool point_is_set(const Mask& mask, int point) {
  return ((mask.words[point / 64] >> (point % 64)) & 1ULL) != 0;
}

void toggle_point(Mask& mask, int point) {
  mask.words[point / 64] ^= 1ULL << (point % 64);
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
  return result;
}

std::uint64_t degree_two_vector(int point) {
  std::array<int, kQuotientDimension> coordinates{};
  for (int variable = 0; variable < kQuotientDimension; ++variable) {
    coordinates[variable] =
        (point >> (kQuotientDimension - 1 - variable)) & 1;
  }
  std::uint64_t result = 1;
  int bit = 0;
  for (const int coordinate : coordinates) {
    ++bit;
    if (coordinate != 0) result |= 1ULL << bit;
  }
  for (int left = 0; left < kQuotientDimension; ++left) {
    for (int right = left + 1; right < kQuotientDimension; ++right) {
      ++bit;
      if (coordinates[left] != 0 && coordinates[right] != 0) {
        result |= 1ULL << bit;
      }
    }
  }
  return result;
}

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 31; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool mask_equal(const Mask& left, const Mask& right) {
  return left.words == right.words;
}

Mask toggle_for_coefficients(std::uint64_t coefficients,
                             const std::vector<int>& core_points) {
  Mask result;
  while (coefficients != 0) {
    const std::uint64_t bit = coefficients & (~coefficients + 1);
    const int point = core_points[__builtin_ctzll(bit)];
    toggle_point(result, 2 * point);
    toggle_point(result, 2 * point + 1);
    coefficients ^= bit;
  }
  return result;
}

void xor_mask(Mask& destination, const Mask& toggle) {
  for (int index = 0; index < 32; ++index) {
    destination.words[index] ^= toggle.words[index];
  }
}

std::array<std::uint64_t, 4> signature(const Mask& mask) {
  std::array<std::uint16_t, kSupportSize> points{};
  int point_count = 0;
  std::array<std::int16_t, kAmbientSize> spectrum{};
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!point_is_set(mask, point)) continue;
    if (point_count >= kSupportSize) {
    throw std::runtime_error("Candidate has weight greater than the configured length");
    }
    points[point_count++] = static_cast<std::uint16_t>(point);
    spectrum[point] = 1;
  }
  if (point_count != kSupportSize) {
    throw std::runtime_error("Candidate does not have the configured weight");
  }

  std::array<std::uint16_t, kAmbientSize> differences{};
  for (int left = 0; left < kSupportSize; ++left) {
    for (int right = left + 1; right < kSupportSize; ++right) {
      differences[points[left] ^ points[right]] += 2;
    }
  }
  for (int step = 1; step < kAmbientSize; step *= 2) {
    for (int block = 0; block < kAmbientSize; block += 2 * step) {
      for (int offset = 0; offset < step; ++offset) {
        const auto left = spectrum[block + offset];
        const auto right = spectrum[block + step + offset];
        spectrum[block + offset] = left + right;
        spectrum[block + step + offset] = left - right;
      }
    }
  }

  std::array<std::uint16_t, kProfileBins> hyperplanes{};
  for (int functional = 1; functional < kAmbientSize; ++functional) {
    const int coefficient = spectrum[functional];
    if ((kSupportSize + coefficient) % 2 != 0) {
      throw std::runtime_error("Impossible Walsh coefficient parity");
    }
    const int intersection = (kSupportSize + coefficient) / 2;
    ++hyperplanes[std::min(intersection, kSupportSize - intersection)];
  }
  std::array<std::uint16_t, kProfileBins> difference_histogram{};
  for (int translation = 1; translation < kAmbientSize; ++translation) {
    const int count = differences[translation];
    if (count % 2 != 0 || count > kSupportSize) {
      throw std::runtime_error("Impossible difference count");
    }
    ++difference_histogram[count / 2];
  }
  std::array<std::array<std::uint8_t, kProfileBins>, kSupportSize> local{};
  for (int row = 0; row < kSupportSize; ++row) {
    for (int column = 0; column < kSupportSize; ++column) {
      if (row != column) {
        ++local[row][differences[points[row] ^ points[column]] / 2];
      }
    }
  }
  std::sort(local.begin(), local.end());

  ProfileHash hash;
  hash.add_byte(0xa0);
  hash.add_u16(kAmbientSize);
  hash.add_u16(kSupportSize);
  hash.add_byte(0xa1);
  for (const auto count : hyperplanes) hash.add_u16(count);
  hash.add_byte(0xa2);
  for (const auto count : difference_histogram) hash.add_u16(count);
  hash.add_byte(0xa3);
  for (const auto& row : local) {
    for (const auto count : row) hash.add_byte(count);
  }
  return hash.finish();
}

std::vector<SourceRecord> read_sources(const Arguments& arguments) {
  std::ifstream input(arguments.input, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open source ledger");
  const std::uint64_t bytes = fs::file_size(arguments.input);
  if (bytes % sizeof(SourceRecord) != 0 ||
      arguments.source_start + arguments.source_count >
          bytes / sizeof(SourceRecord)) {
    throw std::runtime_error("Requested source interval lies outside the ledger");
  }
  input.seekg(static_cast<std::streamoff>(arguments.source_start * sizeof(SourceRecord)));
  std::vector<SourceRecord> records(arguments.source_count);
  input.read(reinterpret_cast<char*>(records.data()),
             static_cast<std::streamsize>(records.size() * sizeof(SourceRecord)));
  if (!input) throw std::runtime_error("Could not read the complete source interval");
  return records;
}

void generate_source(const SourceRecord& source, std::vector<Mask>& candidates) {
  if (source.source_kind != 0 || source.multiplicity != 0) {
    throw std::runtime_error("The m=11 cover only supports N=0 sources");
  }
  const int expected_weight = kSupportSize;
  std::vector<int> core_points;
  std::vector<std::uint64_t> columns;
  core_points.reserve(expected_weight);
  columns.reserve(expected_weight);
  ColumnSolver solver;
  for (int point = 0; point < kQuotientSize; ++point) {
    if (!point_is_set(source.words, point)) continue;
    const int index = static_cast<int>(core_points.size());
    core_points.push_back(point);
    columns.push_back(degree_two_vector(point));
    solver.add(columns.back(), index);
  }
  if (static_cast<int>(core_points.size()) != expected_weight ||
      solver.rank() != source.quadratic_rank) {
    throw std::runtime_error("Source rank or weight mismatch");
  }

  CoefficientBasis quotient_basis;
  std::vector<std::uint64_t> affine_rows;
  affine_rows.push_back((1ULL << expected_weight) - 1);
  for (int coordinate = 0; coordinate < kQuotientDimension; ++coordinate) {
    std::uint64_t row = 0;
    for (int index = 0; index < expected_weight; ++index) {
      if (((core_points[index] >> coordinate) & 1) != 0) row |= 1ULL << index;
    }
    affine_rows.push_back(row);
  }
  for (const auto row : affine_rows) {
    if (!quotient_basis.add(row) || solver.evaluate(row, columns) != 0) {
      throw std::runtime_error("Invalid affine shear row");
    }
  }
  std::vector<std::uint64_t> lift_complement;
  for (const auto relation : solver.relations()) {
    if (quotient_basis.add(relation)) lift_complement.push_back(relation);
  }
  if (quotient_basis.rank() !=
          static_cast<int>(affine_rows.size() + lift_complement.size()) ||
      lift_complement.size() != source.lift_dimension) {
    throw std::runtime_error("Lift-complement profile mismatch");
  }

  Mask singleton_base;
  for (const int point : core_points) toggle_point(singleton_base, 2 * point);
  std::vector<Mask> complement_toggles;
  for (const auto relation : lift_complement) {
    complement_toggles.push_back(toggle_for_coefficients(relation, core_points));
  }
  std::vector<int> fibre_points;
  fibre_points.push_back(-1);
  if (fibre_points.size() != source.compatible_fibre_count) {
    throw std::runtime_error("Compatible-fibre count mismatch");
  }

  const std::uint64_t lift_count = 1ULL << lift_complement.size();
  const std::uint64_t expected_raw = fibre_points.size() * lift_count;
  if (expected_raw != source.raw_marked_extension_count) {
    throw std::runtime_error("Raw extension count mismatch");
  }
  const std::size_t before = candidates.size();
  for (const int fibre_point : fibre_points) {
    std::uint64_t particular = 0;
    const std::uint64_t right_hand_side =
        fibre_point < 0 ? 0 : degree_two_vector(fibre_point);
    if (!solver.solve(right_hand_side, particular)) {
      throw std::runtime_error("A compatible fibre has no particular lift");
    }
    Mask candidate = singleton_base;
    xor_mask(candidate, toggle_for_coefficients(particular, core_points));
    if (fibre_point >= 0) {
      toggle_point(candidate, 2 * fibre_point);
      toggle_point(candidate, 2 * fibre_point + 1);
    }
    for (std::uint64_t coefficient = 0; coefficient < lift_count; ++coefficient) {
      if (source.multiplicity != 0 || coefficient != 0) {
        if (mask_weight(candidate) != kSupportSize) {
          throw std::runtime_error("Generated candidate has incorrect weight");
        }
        candidates.push_back(candidate);
      }
      const std::uint64_t changed = coefficient + 1;
      if (changed < lift_count) {
        xor_mask(candidate, complement_toggles[__builtin_ctzll(changed)]);
      }
    }
  }
  const std::uint64_t expected_full_rank =
      expected_raw - (source.multiplicity == 0 ? 1 : 0);
  if (candidates.size() - before != expected_full_rank) {
    throw std::runtime_error("Full-rank candidate count mismatch");
  }
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input") {
      arguments.input = value;
    } else if (option == "--masks-output") {
      arguments.masks_output = value;
    } else if (option == "--output") {
      arguments.output = value;
    } else if (option == "--source-start") {
      arguments.source_start = std::stoull(value);
    } else if (option == "--source-count") {
      arguments.source_count = std::stoull(value);
    } else if (option == "--threads") {
      arguments.threads = std::stoi(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input.empty() || arguments.masks_output.empty() ||
      arguments.output.empty() ||
      arguments.source_count == 0 || arguments.threads <= 0) {
    throw std::runtime_error(
        "Usage: --input PATH --masks-output PATH --output PATH "
        "--source-start I --source-count N --threads T");
  }
  return arguments;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto sources = read_sources(arguments);
    std::uint64_t expected_occurrences = 0;
    for (const auto& source : sources) {
      expected_occurrences += source.raw_marked_extension_count -
                              (source.multiplicity == 0 ? 1 : 0);
    }
    if (expected_occurrences > std::numeric_limits<std::size_t>::max()) {
      throw std::runtime_error("Candidate count exceeds addressable memory");
    }
    std::vector<Mask> masks;
    masks.reserve(static_cast<std::size_t>(expected_occurrences));
    for (const auto& source : sources) generate_source(source, masks);
    const auto generated = std::chrono::steady_clock::now();
    if (masks.size() != expected_occurrences) {
      throw std::runtime_error("Generated occurrence count mismatch");
    }
    std::sort(masks.begin(), masks.end(), mask_less);
    const auto unique_end = std::unique(masks.begin(), masks.end(), mask_equal);
    const std::uint64_t duplicate_count = masks.end() - unique_end;
    masks.erase(unique_end, masks.end());
    const auto deduplicated = std::chrono::steady_clock::now();

    std::ofstream masks_output(
        arguments.masks_output, std::ios::binary | std::ios::trunc);
    masks_output.write(
        reinterpret_cast<const char*>(masks.data()),
        static_cast<std::streamsize>(masks.size() * sizeof(Mask)));
    if (!masks_output) throw std::runtime_error("Could not write candidate masks");
    masks_output.close();

    std::vector<SignatureRecord> records(masks.size());
#ifdef _OPENMP
    omp_set_num_threads(arguments.threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(masks.size()); ++index) {
      records[index] = SignatureRecord{signature(masks[index]), masks[index]};
    }
    const auto signed_at = std::chrono::steady_clock::now();
    std::sort(records.begin(), records.end(),
              [](const SignatureRecord& left, const SignatureRecord& right) {
                if (left.signature != right.signature) {
                  return left.signature < right.signature;
                }
                return mask_less(left.mask, right.mask);
              });
    std::ofstream output(arguments.output, std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(records.data()),
                 static_cast<std::streamsize>(records.size() * sizeof(SignatureRecord)));
    if (!output) throw std::runtime_error("Could not write signature records");
    output.close();
    const auto finished = std::chrono::steady_clock::now();
    const auto seconds = [](auto first, auto last) {
      return std::chrono::duration<double>(last - first).count();
    };
    std::cout << "{\n"
              << "  \"status\": \"complete_configured_m11_native_fused_signature_batch_v1\",\n"
              << "  \"source_count\": " << sources.size() << ",\n"
              << "  \"full_rank_occurrence_count\": " << expected_occurrences << ",\n"
              << "  \"distinct_candidate_count\": " << records.size() << ",\n"
              << "  \"within_batch_duplicate_count\": " << duplicate_count << ",\n"
              << "  \"threads\": " << arguments.threads << ",\n"
              << "  \"generation_seconds\": " << seconds(started, generated) << ",\n"
              << "  \"deduplication_seconds\": " << seconds(generated, deduplicated)
              << ",\n"
              << "  \"signature_seconds\": " << seconds(deduplicated, signed_at)
              << ",\n"
              << "  \"sort_and_write_seconds\": " << seconds(signed_at, finished)
              << ",\n"
              << "  \"elapsed_seconds\": " << seconds(started, finished) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
