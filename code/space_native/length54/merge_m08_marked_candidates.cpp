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

namespace fs = std::filesystem;

namespace {

struct InputRecord {
  std::array<std::uint64_t, 4> words{};
  std::uint64_t raw_marked_pair_count = 0;
};

struct UnionRecord {
  std::array<std::uint64_t, 4> words{};
  std::uint64_t direction_marked_occurrence_count = 0;
  std::uint64_t raw_marked_pair_count = 0;
};

static_assert(sizeof(InputRecord) == 40);
static_assert(sizeof(UnionRecord) == 48);

struct Arguments {
  fs::path input_directory;
  fs::path output;
  std::uint64_t expected_records = 0;
  std::uint64_t expected_raw_pairs = 0;
  std::uint64_t expected_files = 0;
};

bool mask_less(const InputRecord& left, const InputRecord& right) {
  for (int index = 3; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool same_mask(const InputRecord& left, const InputRecord& right) {
  return left.words == right.words;
}

int support_weight(const std::array<std::uint64_t, 4>& words) {
  int result = 0;
  for (const auto word : words) result += std::popcount(word);
  return result;
}

std::vector<int> support_points(const std::array<std::uint64_t, 4>& words) {
  std::vector<int> points;
  points.reserve(50);
  for (int block = 0; block < 4; ++block) {
    std::uint64_t word = words[block];
    while (word != 0) {
      const int bit = std::countr_zero(word);
      points.push_back(64 * block + bit);
      word &= word - 1;
    }
  }
  return points;
}

int binary_rank(std::vector<std::uint8_t> values) {
  int rank = 0;
  for (int bit = 7; bit >= 0; --bit) {
    auto pivot = std::find_if(
        values.begin() + rank, values.end(),
        [bit](std::uint8_t value) { return ((value >> bit) & 1U) != 0; });
    if (pivot == values.end()) continue;
    std::iter_swap(values.begin() + rank, pivot);
    for (std::size_t row = 0; row < values.size(); ++row) {
      if (static_cast<int>(row) != rank && ((values[row] >> bit) & 1U) != 0) {
        values[row] ^= values[rank];
      }
    }
    ++rank;
  }
  return rank;
}

bool has_full_affine_rank(const std::vector<int>& points) {
  if (points.empty()) return false;
  std::vector<std::uint8_t> differences;
  differences.reserve(points.size() - 1);
  for (std::size_t index = 1; index < points.size(); ++index) {
    differences.push_back(
        static_cast<std::uint8_t>(points[index] ^ points.front()));
  }
  return binary_rank(std::move(differences)) == 8;
}

bool has_even_degree_three_moments(const std::vector<int>& points) {
  std::array<std::uint8_t, 8> singles{};
  std::array<std::array<std::uint8_t, 8>, 8> pairs{};
  std::array<std::array<std::array<std::uint8_t, 8>, 8>, 8> triples{};
  std::uint8_t constant = 0;
  for (const int point : points) {
    constant ^= 1;
    for (int first = 0; first < 8; ++first) {
      if (((point >> first) & 1) == 0) continue;
      singles[first] ^= 1;
      for (int second = first + 1; second < 8; ++second) {
        if (((point >> second) & 1) == 0) continue;
        pairs[first][second] ^= 1;
        for (int third = second + 1; third < 8; ++third) {
          if (((point >> third) & 1) != 0) triples[first][second][third] ^= 1;
        }
      }
    }
  }
  if (constant != 0) return false;
  for (int first = 0; first < 8; ++first) {
    if (singles[first] != 0) return false;
    for (int second = first + 1; second < 8; ++second) {
      if (pairs[first][second] != 0) return false;
      for (int third = second + 1; third < 8; ++third) {
        if (triples[first][second][third] != 0) return false;
      }
    }
  }
  return true;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input-directory") {
      arguments.input_directory = value;
    } else if (option == "--output") {
      arguments.output = value;
    } else if (option == "--expected-records") {
      arguments.expected_records = std::stoull(value);
    } else if (option == "--expected-raw-pairs") {
      arguments.expected_raw_pairs = std::stoull(value);
    } else if (option == "--expected-files") {
      arguments.expected_files = std::stoull(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input_directory.empty() || arguments.output.empty() ||
      arguments.expected_records == 0 || arguments.expected_raw_pairs == 0 ||
      arguments.expected_files == 0) {
    throw std::runtime_error(
        "Required: --input-directory PATH --output PATH --expected-records N "
        "--expected-raw-pairs N --expected-files N");
  }
  return arguments;
}

double elapsed_seconds(
    std::chrono::steady_clock::time_point begin,
    std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_directory(arguments.input_directory)) {
      throw std::runtime_error("Input directory does not exist");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite candidate union");
    }

    std::vector<fs::path> inputs;
    for (const auto& entry : fs::directory_iterator(arguments.input_directory)) {
      if (entry.is_regular_file() && entry.path().extension() == ".bin") {
        inputs.push_back(entry.path());
      }
    }
    std::sort(inputs.begin(), inputs.end());
    if (inputs.size() != arguments.expected_files) {
      throw std::runtime_error("Unexpected marked-parent binary count");
    }

    std::vector<InputRecord> records;
    records.reserve(static_cast<std::size_t>(arguments.expected_records));
    std::uint64_t input_raw_pairs = 0;
    for (const auto& path : inputs) {
      const auto byte_count = fs::file_size(path);
      if (byte_count % sizeof(InputRecord) != 0) {
        throw std::runtime_error("Malformed marked-parent binary: " + path.string());
      }
      const auto old_size = records.size();
      const auto count = byte_count / sizeof(InputRecord);
      records.resize(old_size + static_cast<std::size_t>(count));
      std::ifstream input(path, std::ios::binary);
      input.read(
          reinterpret_cast<char*>(records.data() + old_size),
          static_cast<std::streamsize>(byte_count));
      if (!input || input.gcount() != static_cast<std::streamsize>(byte_count)) {
        throw std::runtime_error("Could not read marked-parent binary: " + path.string());
      }
      for (std::size_t index = old_size; index < records.size(); ++index) {
        if (records[index].raw_marked_pair_count == 0) {
          throw std::runtime_error("A marked orbit has zero multiplicity");
        }
        if (support_weight(records[index].words) != 50) {
          throw std::runtime_error("A marked candidate does not have weight 50");
        }
        input_raw_pairs += records[index].raw_marked_pair_count;
      }
    }
    if (records.size() != arguments.expected_records ||
        input_raw_pairs != arguments.expected_raw_pairs) {
      throw std::runtime_error("Marked-parent aggregate count mismatch");
    }
    const auto read_finished = std::chrono::steady_clock::now();

    std::sort(records.begin(), records.end(), mask_less);
    const auto sort_finished = std::chrono::steady_clock::now();

    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot create candidate union");
    std::uint64_t unique_count = 0;
    std::uint64_t occurrence_count = 0;
    std::uint64_t raw_pair_count = 0;
    std::uint64_t invalid_moment_count = 0;
    std::uint64_t invalid_affine_rank_count = 0;
    for (std::size_t begin = 0; begin < records.size();) {
      std::size_t end = begin + 1;
      std::uint64_t raw = records[begin].raw_marked_pair_count;
      while (end < records.size() && same_mask(records[begin], records[end])) {
        raw += records[end].raw_marked_pair_count;
        ++end;
      }
      const auto points = support_points(records[begin].words);
      if (!has_even_degree_three_moments(points)) ++invalid_moment_count;
      if (!has_full_affine_rank(points)) ++invalid_affine_rank_count;
      const UnionRecord result{
          records[begin].words,
          static_cast<std::uint64_t>(end - begin),
          raw,
      };
      output.write(reinterpret_cast<const char*>(&result), sizeof(result));
      if (!output) throw std::runtime_error("Could not write candidate union");
      ++unique_count;
      occurrence_count += result.direction_marked_occurrence_count;
      raw_pair_count += result.raw_marked_pair_count;
      begin = end;
    }
    output.close();
    if (occurrence_count != arguments.expected_records ||
        raw_pair_count != arguments.expected_raw_pairs ||
        invalid_moment_count != 0 || invalid_affine_rank_count != 0) {
      fs::remove(temporary);
      throw std::runtime_error("Candidate-union validation failed");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();

    std::cout << "{\n"
              << "  \"status\": \"complete_length50_m08_native_candidate_union\",\n"
              << "  \"input_file_count\": " << inputs.size() << ",\n"
              << "  \"direction_marked_occurrence_count\": " << occurrence_count << ",\n"
              << "  \"raw_marked_pair_count\": " << raw_pair_count << ",\n"
              << "  \"distinct_candidate_count\": " << unique_count << ",\n"
              << "  \"invalid_moment_count\": " << invalid_moment_count << ",\n"
              << "  \"invalid_affine_rank_count\": " << invalid_affine_rank_count << ",\n"
              << "  \"read_and_validate_seconds\": "
              << elapsed_seconds(started, read_finished) << ",\n"
              << "  \"sort_seconds\": "
              << elapsed_seconds(read_finished, sort_finished) << ",\n"
              << "  \"aggregate_and_validate_seconds\": "
              << elapsed_seconds(sort_finished, finished) << ",\n"
              << "  \"elapsed_seconds\": " << elapsed_seconds(started, finished) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
