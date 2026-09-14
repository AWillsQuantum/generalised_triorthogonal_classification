#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kQuotientDimension = 8;
constexpr int kTargetDimension = 9;
constexpr int kQuotientSize = 1 << kQuotientDimension;
constexpr int kAmbientSize = 1 << kTargetDimension;
constexpr int kSupportSize = 54;

struct SourceRecord {
  std::array<std::uint64_t, 4> words{};
  std::uint8_t source_kind = 0;
  std::uint8_t multiplicity = 0;
  std::uint8_t quadratic_rank = 0;
  std::uint8_t lift_dimension = 0;
  std::uint32_t reserved = 0;
  std::uint64_t source_index = 0;
  std::uint64_t compatible_fibre_count = 0;
  std::uint64_t raw_marked_extension_count = 0;
};

struct Mask {
  std::array<std::uint64_t, 8> words{};
};

static_assert(sizeof(SourceRecord) == 64);
static_assert(sizeof(Mask) == 64);

struct Arguments {
  fs::path input;
  fs::path output;
  fs::path debug_all_output;
  std::uint64_t source_start = 0;
  std::uint64_t source_count = 0;
  std::uint64_t candidate_limit = 0;
  int expected_multiplicity = -1;
  bool count_only = false;
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

  std::uint64_t evaluate(
      std::uint64_t coefficients,
      const std::vector<std::uint64_t>& columns) const {
    std::uint64_t result = 0;
    while (coefficients != 0) {
      const int bit = __builtin_ctzll(coefficients);
      result ^= columns[bit];
      coefficients &= coefficients - 1;
    }
    return result;
  }

  int rank() const { return rank_; }
  const std::vector<std::uint64_t>& relations() const { return relations_; }

 private:
  void reduce_with_coefficients(
      std::uint64_t& value, std::uint64_t& coefficients) const {
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

struct Accounting {
  std::uint64_t sources_started = 0;
  std::uint64_t sources_completed = 0;
  std::uint64_t expected_raw = 0;
  std::uint64_t expected_full_rank = 0;
  std::uint64_t processed = 0;
  std::array<std::uint64_t, 3> minimum_counts{};
  std::uint64_t retained_occurrences = 0;
};

bool source_point_is_set(
    const std::array<std::uint64_t, 4>& words, int point) {
  return ((words[point / 64] >> (point % 64)) & 1ULL) != 0;
}

void toggle_point(Mask& mask, int point) {
  mask.words[point / 64] ^= 1ULL << (point % 64);
}

void xor_mask(Mask& destination, const Mask& source) {
  for (std::size_t index = 0; index < destination.words.size(); ++index) {
    destination.words[index] ^= source.words[index];
  }
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
  return result;
}

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool mask_equal(const Mask& left, const Mask& right) {
  return left.words == right.words;
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

Mask toggle_for_coefficients(
    std::uint64_t coefficients, const std::vector<int>& core_points) {
  Mask result;
  while (coefficients != 0) {
    const int bit = __builtin_ctzll(coefficients);
    const int point = core_points[bit];
    toggle_point(result, 2 * point);
    toggle_point(result, 2 * point + 1);
    coefficients &= coefficients - 1;
  }
  return result;
}

int minimum_pair_multiplicity(const Mask& mask, int chosen_multiplicity) {
  if (chosen_multiplicity == 0) return 0;
  std::array<std::uint16_t, kSupportSize> points{};
  int point_count = 0;
  for (int word_index = 0; word_index < 8; ++word_index) {
    std::uint64_t word = mask.words[word_index];
    while (word != 0) {
      const int bit = __builtin_ctzll(word);
      if (point_count >= kSupportSize) {
        throw std::runtime_error("Generated candidate exceeds target weight");
      }
      points[point_count++] = static_cast<std::uint16_t>(64 * word_index + bit);
      word &= word - 1;
    }
  }
  if (point_count != kSupportSize) {
    throw std::runtime_error("Generated candidate has incorrect target weight");
  }

  std::array<std::uint64_t, 8> once{};
  std::array<std::uint64_t, 8> twice{};
  for (int left = 0; left < kSupportSize - 1; ++left) {
    for (int right = left + 1; right < kSupportSize; ++right) {
      const int difference = points[left] ^ points[right];
      const int word_index = difference / 64;
      const std::uint64_t bit = 1ULL << (difference % 64);
      twice[word_index] |= once[word_index] & bit;
      once[word_index] |= bit;
    }
  }
  once[0] |= 1ULL;
  twice[0] |= 1ULL;
  const bool has_zero = std::any_of(
      once.begin(), once.end(),
      [](std::uint64_t word) { return word != std::numeric_limits<std::uint64_t>::max(); });
  if (has_zero) return 0;
  if (chosen_multiplicity == 1) return 1;
  const bool has_one = std::any_of(
      twice.begin(), twice.end(),
      [](std::uint64_t word) { return word != std::numeric_limits<std::uint64_t>::max(); });
  return has_one ? 1 : 2;
}

std::vector<SourceRecord> read_sources(const Arguments& arguments) {
  std::ifstream input(arguments.input, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot open the source ledger");
  const std::uint64_t bytes = fs::file_size(arguments.input);
  if (bytes % sizeof(SourceRecord) != 0 ||
      arguments.source_start + arguments.source_count >
          bytes / sizeof(SourceRecord)) {
    throw std::runtime_error("Requested source interval lies outside the ledger");
  }
  input.seekg(
      static_cast<std::streamoff>(arguments.source_start * sizeof(SourceRecord)));
  std::vector<SourceRecord> sources(arguments.source_count);
  input.read(
      reinterpret_cast<char*>(sources.data()),
      static_cast<std::streamsize>(sources.size() * sizeof(SourceRecord)));
  if (!input) throw std::runtime_error("Could not read the complete source interval");
  return sources;
}

bool candidate_limit_reached(
    const Arguments& arguments, const Accounting& accounting) {
  return arguments.candidate_limit != 0 &&
         accounting.processed >= arguments.candidate_limit;
}

bool process_candidate(
    const Mask& candidate, int multiplicity, const Arguments& arguments,
    Accounting& accounting, std::vector<Mask>& retained,
    std::vector<Mask>& debug_all) {
  if (candidate_limit_reached(arguments, accounting)) return false;
  if (mask_weight(candidate) != kSupportSize) {
    throw std::runtime_error("Generated candidate has incorrect target weight");
  }
  const int minimum = minimum_pair_multiplicity(candidate, multiplicity);
  if (minimum < 0 || minimum > multiplicity) {
    throw std::runtime_error("Minimum multiplicity exceeds the chosen direction");
  }
  ++accounting.processed;
  ++accounting.minimum_counts[minimum];
  if (!arguments.debug_all_output.empty()) debug_all.push_back(candidate);
  if (minimum == multiplicity) {
    ++accounting.retained_occurrences;
    if (!arguments.count_only) retained.push_back(candidate);
  }
  return !candidate_limit_reached(arguments, accounting);
}

bool process_fibre(
    const std::vector<int>& fibre_points, const ColumnSolver& solver,
    const std::vector<std::uint64_t>& quadratic,
    const std::vector<int>& core_points, const Mask& singleton_base,
    const std::vector<Mask>& complement_toggles, int multiplicity,
    const Arguments& arguments, Accounting& accounting,
    std::vector<Mask>& retained, std::vector<Mask>& debug_all) {
  std::uint64_t right_hand_side = 0;
  for (const int point : fibre_points) right_hand_side ^= quadratic[point];
  std::uint64_t particular = 0;
  if (!solver.solve(right_hand_side, particular)) {
    throw std::runtime_error("A compatible fibre has no particular lift");
  }
  Mask candidate = singleton_base;
  xor_mask(candidate, toggle_for_coefficients(particular, core_points));
  for (const int point : fibre_points) {
    toggle_point(candidate, 2 * point);
    toggle_point(candidate, 2 * point + 1);
  }

  const std::uint64_t lift_count = 1ULL << complement_toggles.size();
  for (std::uint64_t coefficient = 0; coefficient < lift_count; ++coefficient) {
    if (multiplicity != 0 || coefficient != 0) {
      if (!process_candidate(
              candidate, multiplicity, arguments, accounting, retained,
              debug_all)) {
        return false;
      }
    }
    const std::uint64_t changed = coefficient + 1;
    if (changed < lift_count) {
      xor_mask(candidate, complement_toggles[__builtin_ctzll(changed)]);
    }
  }
  return true;
}

bool generate_source(
    const SourceRecord& source, const Arguments& arguments,
    Accounting& accounting, std::vector<Mask>& retained,
    std::vector<Mask>& debug_all) {
  const int multiplicity = source.multiplicity;
  if (source.source_kind != source.multiplicity || multiplicity < 0 ||
      multiplicity > 2 || source.reserved != 0 ||
      (arguments.expected_multiplicity >= 0 &&
       multiplicity != arguments.expected_multiplicity)) {
    throw std::runtime_error("Source kind, multiplicity, or reserved field is invalid");
  }
  const int expected_weight = kSupportSize - 2 * multiplicity;
  std::vector<int> core_points;
  std::vector<std::uint64_t> columns;
  core_points.reserve(expected_weight);
  columns.reserve(expected_weight);
  ColumnSolver solver;
  for (int point = 0; point < kQuotientSize; ++point) {
    if (!source_point_is_set(source.words, point)) continue;
    const int index = static_cast<int>(core_points.size());
    core_points.push_back(point);
    columns.push_back(degree_two_vector(point));
    solver.add(columns.back(), index);
  }
  if (static_cast<int>(core_points.size()) != expected_weight ||
      solver.rank() != source.quadratic_rank) {
    throw std::runtime_error("Source weight or quadratic rank is invalid");
  }

  CoefficientBasis quotient_basis;
  std::vector<std::uint64_t> affine_rows;
  affine_rows.push_back((1ULL << expected_weight) - 1);
  for (int coordinate = 0; coordinate < kQuotientDimension; ++coordinate) {
    std::uint64_t row = 0;
    for (int index = 0; index < expected_weight; ++index) {
      if (((core_points[index] >> coordinate) & 1) != 0) {
        row |= 1ULL << index;
      }
    }
    affine_rows.push_back(row);
  }
  for (const auto row : affine_rows) {
    if (!quotient_basis.add(row) || solver.evaluate(row, columns) != 0) {
      throw std::runtime_error("Source affine rows are invalid");
    }
  }
  std::vector<std::uint64_t> lift_complement;
  for (const auto relation : solver.relations()) {
    if (quotient_basis.add(relation)) lift_complement.push_back(relation);
  }
  if (quotient_basis.rank() !=
          static_cast<int>(affine_rows.size() + lift_complement.size()) ||
      lift_complement.size() != source.lift_dimension) {
    throw std::runtime_error("Source lift-complement profile is invalid");
  }

  Mask singleton_base;
  for (const int point : core_points) toggle_point(singleton_base, 2 * point);
  std::vector<Mask> complement_toggles;
  complement_toggles.reserve(lift_complement.size());
  for (const auto relation : lift_complement) {
    complement_toggles.push_back(
        toggle_for_coefficients(relation, core_points));
  }

  std::vector<std::uint64_t> quadratic(kQuotientSize);
  std::map<std::uint64_t, std::vector<int>> external_by_label;
  for (int point = 0; point < kQuotientSize; ++point) {
    quadratic[point] = degree_two_vector(point);
    if (!source_point_is_set(source.words, point)) {
      external_by_label[solver.reduce(quadratic[point])].push_back(point);
    }
  }

  std::uint64_t compatible = 0;
  for (const auto& [label, points] : external_by_label) {
    if (multiplicity == 1 && label == 0) compatible += points.size();
    if (multiplicity == 2) {
      compatible += static_cast<std::uint64_t>(points.size()) *
                    (points.size() - 1) / 2;
    }
  }
  if (multiplicity == 0) compatible = 1;
  const std::uint64_t lift_count = 1ULL << lift_complement.size();
  if (compatible != source.compatible_fibre_count ||
      compatible * lift_count != source.raw_marked_extension_count) {
    throw std::runtime_error("Source compatible-fibre accounting is invalid");
  }

  const auto before = accounting.processed;
  bool complete = true;
  if (multiplicity == 0) {
    complete = process_fibre(
        {}, solver, quadratic, core_points, singleton_base,
        complement_toggles, multiplicity, arguments, accounting, retained,
        debug_all);
  } else if (multiplicity == 1) {
    const auto found = external_by_label.find(0);
    if (found == external_by_label.end()) {
      throw std::runtime_error("An N=1 source has no compatible fibre");
    }
    for (const int point : found->second) {
      if (!process_fibre(
              {point}, solver, quadratic, core_points, singleton_base,
              complement_toggles, multiplicity, arguments, accounting,
              retained, debug_all)) {
        complete = false;
        break;
      }
    }
  } else {
    for (const auto& [label, points] : external_by_label) {
      (void)label;
      for (std::size_t left = 0; left + 1 < points.size(); ++left) {
        for (std::size_t right = left + 1; right < points.size(); ++right) {
          if (!process_fibre(
                  {points[left], points[right]}, solver, quadratic,
                  core_points, singleton_base, complement_toggles,
                  multiplicity, arguments, accounting, retained,
                  debug_all)) {
            complete = false;
            break;
          }
        }
        if (!complete) break;
      }
      if (!complete) break;
    }
  }

  if (complete) {
    const std::uint64_t expected =
        source.raw_marked_extension_count - (multiplicity == 0 ? 1 : 0);
    if (accounting.processed - before != expected) {
      throw std::runtime_error("Generated source count does not match its ledger");
    }
  }
  return complete;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--count-only") {
      result.count_only = true;
      continue;
    }
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input") {
      result.input = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--debug-all-output") {
      result.debug_all_output = value;
    } else if (option == "--source-start") {
      result.source_start = std::stoull(value);
    } else if (option == "--source-count") {
      result.source_count = std::stoull(value);
    } else if (option == "--candidate-limit") {
      result.candidate_limit = std::stoull(value);
    } else if (option == "--expected-multiplicity") {
      result.expected_multiplicity = std::stoi(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.input.empty() || result.source_count == 0 ||
      result.expected_multiplicity < 0 || result.expected_multiplicity > 2 ||
      (result.count_only && !result.output.empty()) ||
      (!result.count_only && result.output.empty())) {
    throw std::runtime_error(
        "Usage: --input PATH --source-start I --source-count N "
        "--expected-multiplicity 0|1|2 [--candidate-limit N] "
        "[--debug-all-output PATH] (--count-only | --output PATH)");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto sources = read_sources(arguments);
    const auto loaded = std::chrono::steady_clock::now();
    Accounting accounting;
    for (const auto& source : sources) {
      accounting.expected_raw += source.raw_marked_extension_count;
      accounting.expected_full_rank +=
          source.raw_marked_extension_count -
          (source.multiplicity == 0 ? 1 : 0);
    }
    std::vector<Mask> retained;
    std::vector<Mask> debug_all;
    if (!arguments.count_only && arguments.expected_multiplicity == 0) {
      retained.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(
          accounting.expected_full_rank,
          arguments.candidate_limit == 0 ? accounting.expected_full_rank
                                         : arguments.candidate_limit)));
    }
    if (!arguments.debug_all_output.empty()) {
      debug_all.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(
          accounting.expected_full_rank,
          arguments.candidate_limit == 0 ? accounting.expected_full_rank
                                         : arguments.candidate_limit)));
    }
    bool complete = true;
    for (const auto& source : sources) {
      if (candidate_limit_reached(arguments, accounting)) {
        complete = false;
        break;
      }
      ++accounting.sources_started;
      if (!generate_source(
              source, arguments, accounting, retained, debug_all)) {
        complete = false;
        break;
      }
      ++accounting.sources_completed;
    }
    if (complete && accounting.sources_completed != sources.size()) {
      throw std::runtime_error("Source completion accounting is inconsistent");
    }
    if (complete && accounting.processed != accounting.expected_full_rank) {
      throw std::runtime_error("Batch candidate accounting is inconsistent");
    }
    if (accounting.processed !=
        accounting.minimum_counts[0] + accounting.minimum_counts[1] +
            accounting.minimum_counts[2]) {
      throw std::runtime_error("Minimum-multiplicity accounting is inconsistent");
    }
    const auto generated = std::chrono::steady_clock::now();

    const std::uint64_t retained_occurrences = accounting.retained_occurrences;
    std::uint64_t duplicate_count = 0;
    std::uint64_t debug_duplicate_count = 0;
    if (!arguments.count_only) {
      std::sort(retained.begin(), retained.end(), mask_less);
      const auto unique_end = std::unique(retained.begin(), retained.end(), mask_equal);
      duplicate_count = retained.end() - unique_end;
      retained.erase(unique_end, retained.end());
      if (fs::exists(arguments.output)) {
        throw std::runtime_error("Refusing to overwrite a filtered-mask output");
      }
      const auto temporary = arguments.output.string() + ".tmp";
      std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
      output.write(
          reinterpret_cast<const char*>(retained.data()),
          static_cast<std::streamsize>(retained.size() * sizeof(Mask)));
      output.close();
      if (!output || fs::file_size(temporary) != retained.size() * sizeof(Mask)) {
        fs::remove(temporary);
        throw std::runtime_error("Could not write the filtered-mask output");
      }
      fs::rename(temporary, arguments.output);
    }
    if (!arguments.debug_all_output.empty()) {
      std::sort(debug_all.begin(), debug_all.end(), mask_less);
      const auto unique_end =
          std::unique(debug_all.begin(), debug_all.end(), mask_equal);
      debug_duplicate_count = debug_all.end() - unique_end;
      debug_all.erase(unique_end, debug_all.end());
      if (fs::exists(arguments.debug_all_output)) {
        throw std::runtime_error("Refusing to overwrite a debug-mask output");
      }
      const auto temporary = arguments.debug_all_output.string() + ".tmp";
      std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
      output.write(
          reinterpret_cast<const char*>(debug_all.data()),
          static_cast<std::streamsize>(debug_all.size() * sizeof(Mask)));
      output.close();
      if (!output ||
          fs::file_size(temporary) != debug_all.size() * sizeof(Mask)) {
        fs::remove(temporary);
        throw std::runtime_error("Could not write the debug-mask output");
      }
      fs::rename(temporary, arguments.debug_all_output);
    }
    const auto finished = std::chrono::steady_clock::now();
    const auto seconds = [](auto first, auto last) {
      return std::chrono::duration<double>(last - first).count();
    };
    std::cout
        << "{\n"
        << "  \"status\": \""
        << (complete ? "complete" : "partial")
        << "_length54_m09_minimum_filter_batch_v1\",\n"
        << "  \"source_start\": " << arguments.source_start << ",\n"
        << "  \"source_count\": " << arguments.source_count << ",\n"
        << "  \"sources_started\": " << accounting.sources_started << ",\n"
        << "  \"sources_completed\": " << accounting.sources_completed << ",\n"
        << "  \"multiplicity\": " << arguments.expected_multiplicity << ",\n"
        << "  \"expected_raw_marked_extension_count\": "
        << accounting.expected_raw << ",\n"
        << "  \"expected_full_rank_marked_extension_count\": "
        << accounting.expected_full_rank << ",\n"
        << "  \"processed_candidate_count\": " << accounting.processed << ",\n"
        << "  \"minimum_multiplicity_counts\": ["
        << accounting.minimum_counts[0] << ", "
        << accounting.minimum_counts[1] << ", "
        << accounting.minimum_counts[2] << "],\n"
        << "  \"retained_occurrence_count\": " << retained_occurrences << ",\n"
        << "  \"retained_distinct_batch_mask_count\": " << retained.size() << ",\n"
        << "  \"within_batch_duplicate_count\": " << duplicate_count << ",\n"
        << "  \"debug_all_distinct_batch_mask_count\": " << debug_all.size()
        << ",\n"
        << "  \"debug_all_within_batch_duplicate_count\": "
        << debug_duplicate_count << ",\n"
        << "  \"count_only\": " << (arguments.count_only ? "true" : "false")
        << ",\n"
        << "  \"candidate_limit\": " << arguments.candidate_limit << ",\n"
        << "  \"mask_record_format\": \"<8Q\",\n"
        << "  \"mask_record_bytes\": " << sizeof(Mask) << ",\n"
        << "  \"timings_seconds\": {\n"
        << "    \"load\": " << seconds(started, loaded) << ",\n"
        << "    \"generation_and_filter\": " << seconds(loaded, generated) << ",\n"
        << "    \"sort_deduplicate_and_write\": " << seconds(generated, finished)
        << ",\n"
        << "    \"total\": " << seconds(started, finished) << "\n"
        << "  }\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "m09_minimum_filter_kernel: " << error.what() << "\n";
    return 1;
  }
}
