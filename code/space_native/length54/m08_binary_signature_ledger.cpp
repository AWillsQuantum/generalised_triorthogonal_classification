#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr int kAmbientSize = 256;
constexpr int kSupportSize = 54;
constexpr int kProfileBins = kSupportSize / 2 + 1;

struct UnionRecord {
  std::array<std::uint64_t, 4> mask{};
  std::uint64_t direction_marked_occurrence_count = 0;
  std::uint64_t raw_marked_pair_count = 0;
};

struct MarkedRecord {
  std::array<std::uint64_t, 4> mask{};
  std::uint64_t raw_marked_pair_count = 0;
};

struct LedgerRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 4> mask{};
};

static_assert(sizeof(UnionRecord) == 48);
static_assert(sizeof(MarkedRecord) == 40);
static_assert(sizeof(LedgerRecord) == 64);

struct Arguments {
  fs::path input;
  fs::path output;
  std::string input_format;
  std::uint64_t expected_records = 0;
  std::uint64_t expected_occurrences = 0;
  std::uint64_t expected_raw_pairs = 0;
  int threads = 1;
  std::uint64_t offset = 0;
  std::uint64_t limit = 0;
};

class ProfileHash {
 public:
  ProfileHash()
      : state_{0xcbf29ce484222325ULL, 0x84222325cbf29ce4ULL,
               0x9e3779b97f4a7c15ULL, 0x6a09e667f3bcc909ULL} {}

  void add(std::uint8_t value) {
    constexpr std::array<std::uint64_t, 4> primes{
        0x00000100000001b3ULL, 0x9e3779b185ebca87ULL,
        0xc2b2ae3d27d4eb4fULL, 0x165667b19e3779f9ULL};
    for (std::size_t index = 0; index < state_.size(); ++index) {
      state_[index] ^= static_cast<std::uint64_t>(value) +
                       0x9dU * static_cast<std::uint64_t>(index);
      state_[index] *= primes[index];
      state_[index] ^= state_[index] >> (29 + index);
    }
  }

  std::array<std::uint64_t, 4> finish() const { return state_; }

 private:
  std::array<std::uint64_t, 4> state_;
};

bool point_is_set(const std::array<std::uint64_t, 4>& mask, int point) {
  return ((mask[point / 64] >> (point % 64)) & 1U) != 0;
}

std::array<std::uint64_t, 4> signature(
    const std::array<std::uint64_t, 4>& mask) {
  std::array<std::uint8_t, kSupportSize> points{};
  int point_count = 0;
  std::array<std::int16_t, kAmbientSize> spectrum{};
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!point_is_set(mask, point)) continue;
    if (point_count == kSupportSize) {
      throw std::runtime_error("Candidate has weight greater than 54");
    }
    points[point_count++] = static_cast<std::uint8_t>(point);
    spectrum[point] = 1;
  }
  if (point_count != kSupportSize) {
    throw std::runtime_error("Candidate does not have weight 54");
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
        const std::int16_t left = spectrum[block + offset];
        const std::int16_t right = spectrum[block + step + offset];
        spectrum[block + offset] = left + right;
        spectrum[block + step + offset] = left - right;
      }
    }
  }

  std::array<std::uint8_t, kProfileBins> hyperplanes{};
  for (int functional = 1; functional < kAmbientSize; ++functional) {
    const int coefficient = spectrum[functional];
    if ((kSupportSize + coefficient) % 2 != 0) {
      throw std::runtime_error("Impossible Walsh coefficient parity");
    }
    const int intersection = (kSupportSize + coefficient) / 2;
    const int bin = std::min(intersection, kSupportSize - intersection);
    ++hyperplanes[bin];
  }

  std::array<std::uint8_t, kProfileBins> difference_histogram{};
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
      if (row == column) continue;
      const int count = differences[points[row] ^ points[column]];
      ++local[row][count / 2];
    }
  }
  std::sort(local.begin(), local.end());

  ProfileHash hash;
  hash.add(0xa1);
  for (const auto count : hyperplanes) hash.add(count);
  hash.add(0xa2);
  for (const auto count : difference_histogram) hash.add(count);
  hash.add(0xa3);
  for (const auto& row : local) {
    for (const auto count : row) hash.add(count);
  }
  return hash.finish();
}

bool mask_less(
    const std::array<std::uint64_t, 4>& left,
    const std::array<std::uint64_t, 4>& right) {
  for (int index = 3; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool ledger_less(const LedgerRecord& left, const LedgerRecord& right) {
  if (left.signature != right.signature) {
    return left.signature < right.signature;
  }
  return mask_less(left.mask, right.mask);
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input") {
      arguments.input = value;
    } else if (option == "--output") {
      arguments.output = value;
    } else if (option == "--input-format") {
      arguments.input_format = value;
    } else if (option == "--expected-records") {
      arguments.expected_records = std::stoull(value);
    } else if (option == "--expected-occurrences") {
      arguments.expected_occurrences = std::stoull(value);
    } else if (option == "--expected-raw-pairs") {
      arguments.expected_raw_pairs = std::stoull(value);
    } else if (option == "--threads") {
      arguments.threads = std::stoi(value);
    } else if (option == "--offset-records") {
      arguments.offset = std::stoull(value);
    } else if (option == "--limit") {
      arguments.limit = std::stoull(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input.empty() || arguments.output.empty() ||
      (arguments.input_format != "union" &&
       arguments.input_format != "marked") ||
      arguments.expected_records == 0 || arguments.expected_occurrences == 0 ||
      arguments.threads < 1) {
    throw std::runtime_error(
        "Required: --input PATH --input-format union|marked --output PATH "
        "--expected-records N "
        "--expected-occurrences N --expected-raw-pairs N --threads N "
        "[--offset-records N --limit N]");
  }
  return arguments;
}

double seconds(
    std::chrono::steady_clock::time_point begin,
    std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.input)) {
      throw std::runtime_error("Candidate union does not exist");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite signature ledger");
    }
    const auto input_bytes = fs::file_size(arguments.input);
    const std::uint64_t input_record_bytes =
        arguments.input_format == "union" ? sizeof(UnionRecord)
                                           : sizeof(MarkedRecord);
    if (input_bytes % input_record_bytes != 0) {
      throw std::runtime_error("Candidate union has the wrong byte length");
    }
    const auto available_records = input_bytes / input_record_bytes;
    if (arguments.offset > available_records) {
      throw std::runtime_error("Candidate slice offset exceeds the input ledger");
    }
    const auto remaining_records = available_records - arguments.offset;
    const auto record_count = arguments.limit == 0
                                  ? remaining_records
                                  : std::min(arguments.limit, remaining_records);
    if (record_count != arguments.expected_records) {
      throw std::runtime_error("Candidate slice record count mismatch");
    }

    std::vector<LedgerRecord> ledger(static_cast<std::size_t>(record_count));
    std::ifstream input(arguments.input, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(arguments.offset * input_record_bytes));
    if (!input) throw std::runtime_error("Could not seek to candidate slice offset");
    std::uint64_t occurrence_count = 0;
    std::uint64_t raw_pair_count = 0;
    std::array<std::uint64_t, 4> prior_mask{};
    bool have_prior = false;
    for (std::size_t index = 0; index < ledger.size(); ++index) {
      if (arguments.input_format == "union") {
        UnionRecord record;
        input.read(reinterpret_cast<char*>(&record), sizeof(record));
        if (!input) throw std::runtime_error("Candidate union is truncated");
        if (record.direction_marked_occurrence_count == 0 ||
            record.raw_marked_pair_count == 0) {
          throw std::runtime_error("Candidate union has a zero multiplicity");
        }
        if (have_prior && !mask_less(prior_mask, record.mask)) {
          throw std::runtime_error("Candidate union is not strictly mask-sorted");
        }
        have_prior = true;
        prior_mask = record.mask;
        ledger[index].mask = record.mask;
        occurrence_count += record.direction_marked_occurrence_count;
        raw_pair_count += record.raw_marked_pair_count;
      } else {
        MarkedRecord record;
        input.read(reinterpret_cast<char*>(&record), sizeof(record));
        if (!input) throw std::runtime_error("Marked ledger is truncated");
        if (record.raw_marked_pair_count == 0) {
          throw std::runtime_error("Marked ledger has a zero orbit size");
        }
        ledger[index].mask = record.mask;
        ++occurrence_count;
        raw_pair_count += record.raw_marked_pair_count;
      }
    }
    if (occurrence_count != arguments.expected_occurrences ||
        (arguments.expected_raw_pairs != 0 &&
         raw_pair_count != arguments.expected_raw_pairs)) {
      throw std::runtime_error("Candidate union multiplicity totals mismatch");
    }
    const auto read_finished = std::chrono::steady_clock::now();

#ifdef _OPENMP
    omp_set_num_threads(arguments.threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(ledger.size()); ++index) {
      ledger[static_cast<std::size_t>(index)].signature =
          signature(ledger[static_cast<std::size_t>(index)].mask);
    }
    const auto compute_finished = std::chrono::steady_clock::now();
    std::sort(ledger.begin(), ledger.end(), ledger_less);
    const auto sort_finished = std::chrono::steady_clock::now();

    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    output.write(
        reinterpret_cast<const char*>(ledger.data()),
        static_cast<std::streamsize>(ledger.size() * sizeof(LedgerRecord)));
    output.close();
    if (!output || fs::file_size(temporary) != ledger.size() * sizeof(LedgerRecord)) {
      fs::remove(temporary);
      throw std::runtime_error("Could not write signature ledger");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();
    const double compute_seconds = seconds(read_finished, compute_finished);
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m08_binary_signature_ledger_v1\",\n"
              << "  \"candidate_count\": " << ledger.size() << ",\n"
              << "  \"input_offset_records\": " << arguments.offset << ",\n"
              << "  \"input_format\": \"" << arguments.input_format << "\",\n"
              << "  \"direction_marked_occurrence_count\": "
              << occurrence_count << ",\n"
              << "  \"raw_marked_pair_count\": " << raw_pair_count << ",\n"
              << "  \"threads\": " << arguments.threads << ",\n"
              << "  \"read_seconds\": " << seconds(started, read_finished) << ",\n"
              << "  \"compute_seconds\": " << compute_seconds << ",\n"
              << "  \"sort_seconds\": " << seconds(compute_finished, sort_finished) << ",\n"
              << "  \"write_seconds\": " << seconds(sort_finished, finished) << ",\n"
              << "  \"elapsed_seconds\": " << seconds(started, finished) << ",\n"
              << "  \"candidates_per_compute_second\": "
              << (compute_seconds == 0.0 ? 0.0 : ledger.size() / compute_seconds)
              << "\n}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
