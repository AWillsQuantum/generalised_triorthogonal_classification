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

using Signature = std::array<std::uint64_t, 4>;
using Mask = std::array<std::uint64_t, 8>;

struct WeightedRecord {
  Signature signature{};
  Mask mask{};
  std::uint64_t members = 0;
};

struct ProjectedRecord {
  Signature signature{};
  Mask mask{};
};

static_assert(sizeof(WeightedRecord) == 104);
static_assert(sizeof(ProjectedRecord) == 96);

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool less(const WeightedRecord& left, const WeightedRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

template <typename Record>
bool read_record(std::ifstream& input, Record& record, const char* label) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error(std::string(label) + " is truncated");
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    fs::path weighted_path;
    fs::path projected_path;
    for (int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if (++index >= argc) throw std::runtime_error("Missing option value");
      const fs::path value = argv[index];
      if (option == "--weighted-ledger") {
        weighted_path = value;
      } else if (option == "--projected-ledger") {
        projected_path = value;
      } else {
        throw std::runtime_error("Unknown option: " + option);
      }
    }
    if (weighted_path.empty() || projected_path.empty()) {
      throw std::runtime_error(
          "Required: --weighted-ledger PATH --projected-ledger PATH");
    }
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(weighted_path) ||
        fs::file_size(weighted_path) % sizeof(WeightedRecord) ||
        !fs::is_regular_file(projected_path) ||
        fs::file_size(projected_path) % sizeof(ProjectedRecord)) {
      throw std::runtime_error("A projection audit input has the wrong size");
    }
    std::ifstream weighted(weighted_path, std::ios::binary);
    std::ifstream projected(projected_path, std::ios::binary);
    if (!weighted || !projected) throw std::runtime_error("Could not open audit inputs");
    WeightedRecord source{};
    WeightedRecord prior{};
    ProjectedRecord target{};
    bool have_prior = false;
    std::uint64_t records = 0;
    std::uint64_t members = 0;
    while (read_record(weighted, source, "The weighted ledger")) {
      if (!read_record(projected, target, "The projected ledger") ||
          source.members == 0 || (have_prior && !less(prior, source)) ||
          source.signature != target.signature || source.mask != target.mask ||
          std::numeric_limits<std::uint64_t>::max() - members < source.members) {
        throw std::runtime_error("A projected record is incorrect");
      }
      prior = source;
      have_prior = true;
      members += source.members;
      ++records;
    }
    ProjectedRecord trailing;
    if (read_record(projected, trailing, "The projected ledger") ||
        fs::file_size(projected_path) != records * sizeof(ProjectedRecord)) {
      throw std::runtime_error("The projection audit did not close");
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_weighted_ledger_projection_v1\",\n"
        << "  \"record_count\": " << records << ",\n"
        << "  \"source_member_count\": " << members << ",\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(weighted_path) << ",\n"
        << "  \"projected_ledger_bytes\": " << fs::file_size(projected_path) << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
