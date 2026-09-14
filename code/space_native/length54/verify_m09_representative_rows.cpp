#include <array>
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

struct RepresentativeRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 8> mask{};
  std::uint64_t member_count = 0;
};

static_assert(sizeof(RepresentativeRecord) == 104);

int binary_rank(std::vector<std::uint16_t> values) {
  int rank = 0;
  for (int bit = 9; bit >= 0; --bit) {
    std::size_t pivot = rank;
    while (pivot < values.size() && ((values[pivot] >> bit) & 1U) == 0) {
      ++pivot;
    }
    if (pivot == values.size()) continue;
    std::swap(values[rank], values[pivot]);
    for (std::size_t row = 0; row < values.size(); ++row) {
      if (row != static_cast<std::size_t>(rank) &&
          ((values[row] >> bit) & 1U)) {
        values[row] ^= values[rank];
      }
    }
    ++rank;
  }
  return rank;
}

std::array<std::uint64_t, 10> generator_rows(
    const RepresentativeRecord& record, std::vector<std::uint16_t>& columns) {
  columns.clear();
  columns.reserve(54);
  std::array<std::uint64_t, 10> rows{};
  int column = 0;
  for (int point = 0; point < 512; ++point) {
    if (((record.mask[point / 64] >> (point % 64)) & 1U) == 0) continue;
    if (column == 54) throw std::runtime_error("Representative weight exceeds 54");
    columns.push_back(static_cast<std::uint16_t>(1U | (point << 1)));
    rows[0] |= std::uint64_t{1} << column;
    for (int bit = 0; bit < 9; ++bit) {
      if ((point >> bit) & 1U) rows[bit + 1] |= std::uint64_t{1} << column;
    }
    ++column;
  }
  if (column != 54) throw std::runtime_error("Representative weight is below 54");
  return rows;
}

bool has_even_moments(const std::array<std::uint64_t, 10>& rows) {
  for (int left = 0; left < 10; ++left) {
    for (int middle = left; middle < 10; ++middle) {
      for (int right = middle; right < 10; ++right) {
        if (__builtin_popcountll(rows[left] & rows[middle] & rows[right]) & 1U) {
          return false;
        }
      }
    }
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      throw std::runtime_error("Required: REPRESENTATIVE_LEDGER EXPECTED_RECORDS");
    }
    const fs::path input_path = argv[1];
    const std::uint64_t expected_records = std::stoull(argv[2]);
    if (!fs::is_regular_file(input_path) ||
        fs::file_size(input_path) != expected_records * sizeof(RepresentativeRecord)) {
      throw std::runtime_error("Representative ledger size is invalid");
    }
    const auto started = std::chrono::steady_clock::now();
    std::ifstream input(input_path, std::ios::binary);
    std::uint64_t invalid_weight_count = 0;
    std::uint64_t invalid_rank_count = 0;
    std::uint64_t invalid_moment_count = 0;
    std::uint64_t zero_member_count = 0;
    std::uint64_t total_members = 0;
    std::vector<std::uint16_t> columns;
    for (std::uint64_t index = 0; index < expected_records; ++index) {
      RepresentativeRecord record;
      input.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (!input) throw std::runtime_error("Representative ledger is truncated");
      try {
        const auto rows = generator_rows(record, columns);
        invalid_rank_count += binary_rank(columns) != 10;
        invalid_moment_count += !has_even_moments(rows);
      } catch (const std::runtime_error&) {
        ++invalid_weight_count;
      }
      zero_member_count += record.member_count == 0;
      total_members += record.member_count;
    }
    RepresentativeRecord extra;
    if (input.read(reinterpret_cast<char*>(&extra), sizeof(extra))) {
      throw std::runtime_error("Representative ledger has trailing records");
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\n"
              << "  \"status\": \"verified_length54_m09_all_representative_rows_v1\",\n"
              << "  \"representative_count\": " << expected_records << ",\n"
              << "  \"invalid_weight_count\": " << invalid_weight_count << ",\n"
              << "  \"invalid_rank_count\": " << invalid_rank_count << ",\n"
              << "  \"invalid_moment_count\": " << invalid_moment_count << ",\n"
              << "  \"zero_member_count\": " << zero_member_count << ",\n"
              << "  \"source_member_count\": " << total_members << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << ",\n"
              << "  \"representatives_per_second\": "
              << (elapsed == 0 ? 0 : expected_records / elapsed) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
