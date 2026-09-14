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

constexpr int kSupportSize = 54;
constexpr int kProfileBins = kSupportSize / 2 + 1;

struct Arguments {
  fs::path input;
  fs::path output;
  int dimension = 0;
  int threads = 1;
};

struct ProfileHash {
  std::array<std::uint64_t, 4> state{
      0xcbf29ce484222325ULL,
      0x84222325cbf29ce4ULL,
      0x9e3779b97f4a7c15ULL,
      0x6a09e667f3bcc909ULL,
  };
  static constexpr std::array<std::uint64_t, 4> primes{
      0x00000100000001b3ULL,
      0x9e3779b185ebca87ULL,
      0xc2b2ae3d27d4eb4fULL,
      0x165667b19e3779f9ULL,
  };

  void add_byte(std::uint8_t value) {
    for (int index = 0; index < 4; ++index) {
      state[index] ^= static_cast<std::uint64_t>(value) + 0x9dULL * index;
      state[index] *= primes[index];
      state[index] ^= state[index] >> (29 + index);
    }
  }

  void add_u16(std::uint16_t value) {
    add_byte(static_cast<std::uint8_t>(value));
    add_byte(static_cast<std::uint8_t>(value >> 8));
  }

  std::array<std::uint64_t, 4> finish() const { return state; }
};

struct Record {
  std::array<std::uint64_t, 4> signature{};
  std::vector<std::uint64_t> mask;
};

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
    } else if (option == "--dimension") {
      arguments.dimension = std::stoi(value);
    } else if (option == "--threads") {
      arguments.threads = std::stoi(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input.empty() || arguments.output.empty() ||
      arguments.dimension < 6 || arguments.dimension > 14 ||
      arguments.threads < 1) {
    throw std::runtime_error(
        "Usage: --input PATH --output PATH --dimension M --threads T");
  }
  return arguments;
}

bool mask_less(const std::vector<std::uint64_t>& left,
               const std::vector<std::uint64_t>& right) {
  for (std::size_t index = left.size(); index-- > 0;) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

std::array<std::uint64_t, 4> signature(
    const std::vector<std::uint64_t>& mask, int dimension) {
  const int ambient_size = 1 << dimension;
  std::array<std::uint16_t, kSupportSize> points{};
  int point_count = 0;
  std::vector<std::int16_t> spectrum(ambient_size);
  for (int word_index = 0; word_index < static_cast<int>(mask.size()); ++word_index) {
    auto word = mask[word_index];
    while (word != 0) {
      const int bit = __builtin_ctzll(word);
      const int point = 64 * word_index + bit;
      if (point >= ambient_size || point_count >= kSupportSize) {
        throw std::runtime_error("A candidate mask has invalid support");
      }
      points[point_count++] = static_cast<std::uint16_t>(point);
      spectrum[point] = 1;
      word &= word - 1;
    }
  }
  if (point_count != kSupportSize) {
    throw std::runtime_error("A candidate mask does not have weight 54");
  }

  std::vector<std::uint16_t> differences(ambient_size);
  for (int left = 0; left < kSupportSize; ++left) {
    for (int right = left + 1; right < kSupportSize; ++right) {
      differences[points[left] ^ points[right]] += 2;
    }
  }
  for (int step = 1; step < ambient_size; step *= 2) {
    for (int block = 0; block < ambient_size; block += 2 * step) {
      for (int offset = 0; offset < step; ++offset) {
        const auto left = spectrum[block + offset];
        const auto right = spectrum[block + step + offset];
        spectrum[block + offset] = left + right;
        spectrum[block + step + offset] = left - right;
      }
    }
  }

  std::array<std::uint16_t, kProfileBins> hyperplanes{};
  for (int functional = 1; functional < ambient_size; ++functional) {
    const int coefficient = spectrum[functional];
    if ((kSupportSize + coefficient) % 2 != 0) {
      throw std::runtime_error("Impossible Walsh coefficient parity");
    }
    const int intersection = (kSupportSize + coefficient) / 2;
    ++hyperplanes[std::min(intersection, kSupportSize - intersection)];
  }
  std::array<std::uint16_t, kProfileBins> difference_histogram{};
  for (int translation = 1; translation < ambient_size; ++translation) {
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
  hash.add_u16(static_cast<std::uint16_t>(ambient_size));
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

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const std::size_t word_count = (1ULL << arguments.dimension) / 64;
    const std::uint64_t record_bytes = word_count * sizeof(std::uint64_t);
    if (!fs::is_regular_file(arguments.input) || fs::exists(arguments.output) ||
        fs::file_size(arguments.input) % record_bytes != 0) {
      throw std::runtime_error("Signature input or output path is invalid");
    }
    const std::uint64_t record_count = fs::file_size(arguments.input) / record_bytes;
    std::vector<Record> records(record_count);
    std::ifstream input(arguments.input, std::ios::binary);
    for (auto& record : records) {
      record.mask.resize(word_count);
      input.read(reinterpret_cast<char*>(record.mask.data()),
                 static_cast<std::streamsize>(record_bytes));
      if (!input) throw std::runtime_error("Candidate-mask input is truncated");
    }
    const auto loaded = std::chrono::steady_clock::now();
#ifdef _OPENMP
    omp_set_num_threads(arguments.threads);
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index = 0;
         index < static_cast<std::int64_t>(records.size()); ++index) {
      records[index].signature = signature(records[index].mask, arguments.dimension);
    }
    const auto signed_at = std::chrono::steady_clock::now();
    std::sort(records.begin(), records.end(), [](const Record& left, const Record& right) {
      if (left.signature != right.signature) return left.signature < right.signature;
      return mask_less(left.mask, right.mask);
    });
    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    for (const auto& record : records) {
      output.write(reinterpret_cast<const char*>(record.signature.data()),
                   4 * sizeof(std::uint64_t));
      output.write(reinterpret_cast<const char*>(record.mask.data()),
                   static_cast<std::streamsize>(record_bytes));
    }
    output.close();
    if (!output || fs::file_size(temporary) !=
                       record_count * (record_bytes + 4 * sizeof(std::uint64_t))) {
      fs::remove(temporary);
      throw std::runtime_error("Signature output accounting failed");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();
    const auto seconds = [](auto first, auto last) {
      return std::chrono::duration<double>(last - first).count();
    };
    std::cout << "{\n"
              << "  \"status\": \"complete_runtime_affine_signature\",\n"
              << "  \"dimension\": " << arguments.dimension << ",\n"
              << "  \"record_count\": " << record_count << ",\n"
              << "  \"threads\": " << arguments.threads << ",\n"
              << "  \"load_seconds\": " << seconds(started, loaded) << ",\n"
              << "  \"signature_seconds\": " << seconds(loaded, signed_at) << ",\n"
              << "  \"sort_and_write_seconds\": " << seconds(signed_at, finished) << ",\n"
              << "  \"elapsed_seconds\": " << seconds(started, finished) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
