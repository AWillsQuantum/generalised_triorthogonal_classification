#include "configuration.hpp"
#include <algorithm>
#include <array>
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

constexpr int kAmbientSize = 1 << UTSP_DIMENSION;
constexpr int kSupportSize = UTSP_LENGTH;
constexpr int kProfileBins = kSupportSize / 2 + 1;

struct Mask {
  std::array<std::uint64_t, kAmbientSize / 64> words{};
};

struct Record {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

static_assert(sizeof(Mask) == kAmbientSize / 8);
static_assert(sizeof(Record) == 32 + kAmbientSize / 8);

struct Arguments {
  fs::path input;
  fs::path output;
  int threads = 1;
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
      state_[index] ^= static_cast<std::uint64_t>(value) +
                       0x9dU * static_cast<std::uint64_t>(index);
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

bool point_is_set(const Mask& mask, int point) {
  return ((mask.words[point / 64] >> (point % 64)) & 1U) != 0;
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
  return result;
}

std::array<std::uint64_t, 4> signature(const Mask& mask) {
  std::array<std::uint16_t, kSupportSize> points{};
  int point_count = 0;
  std::array<std::int16_t, kAmbientSize> spectrum{};
  for (int point = 0; point < kAmbientSize; ++point) {
    if (!point_is_set(mask, point)) continue;
    if (point_count == kSupportSize) {
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
      if (row == column) continue;
      ++local[row][differences[points[row] ^ points[column]] / 2];
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

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = kAmbientSize / 64 - 1; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool record_less(const Record& left, const Record& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
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
    } else if (option == "--threads") {
      arguments.threads = std::stoi(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input.empty() || arguments.output.empty() || arguments.threads < 1) {
    throw std::runtime_error("Required: --input PATH --output PATH --threads N");
  }
  return arguments;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.input) ||
        fs::file_size(arguments.input) % sizeof(Mask) != 0) {
      throw std::runtime_error("Candidate input has the wrong byte length");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite signature output");
    }
    const auto record_count = fs::file_size(arguments.input) / sizeof(Mask);
    std::vector<Record> records(static_cast<std::size_t>(record_count));
    std::ifstream input(arguments.input, std::ios::binary);
    for (auto& record : records) {
      input.read(reinterpret_cast<char*>(&record.mask), sizeof(record.mask));
      if (!input || mask_weight(record.mask) != kSupportSize) {
        throw std::runtime_error("Malformed candidate input");
      }
    }
    const auto read_finished = std::chrono::steady_clock::now();
#ifdef _OPENMP
    omp_set_num_threads(arguments.threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(records.size()); ++index) {
      records[static_cast<std::size_t>(index)].signature =
          signature(records[static_cast<std::size_t>(index)].mask);
    }
    const auto compute_finished = std::chrono::steady_clock::now();
    std::sort(records.begin(), records.end(), record_less);
    const auto sort_finished = std::chrono::steady_clock::now();
    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    output.write(
        reinterpret_cast<const char*>(records.data()),
        static_cast<std::streamsize>(records.size() * sizeof(Record)));
    output.close();
    if (!output || fs::file_size(temporary) != records.size() * sizeof(Record)) {
      fs::remove(temporary);
      throw std::runtime_error("Could not write signature output");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();
    const auto seconds = [](auto begin, auto end) {
      return std::chrono::duration<double>(end - begin).count();
    };
    const double compute_seconds = seconds(read_finished, compute_finished);
    std::cout << "{\n"
              << "  \"status\": \"complete_configured_space_signature_v1\",\n"
              << "  \"candidate_count\": " << records.size() << ",\n"
              << "  \"threads\": " << arguments.threads << ",\n"
              << "  \"read_seconds\": " << seconds(started, read_finished) << ",\n"
              << "  \"compute_seconds\": " << compute_seconds << ",\n"
              << "  \"sort_seconds\": " << seconds(compute_finished, sort_finished) << ",\n"
              << "  \"write_seconds\": " << seconds(sort_finished, finished) << ",\n"
              << "  \"candidates_per_compute_second\": "
              << (compute_seconds == 0 ? 0 : records.size() / compute_seconds) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
