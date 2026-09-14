#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>

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

}  // namespace

int main(int argc, char** argv) {
  fs::path input_path;
  fs::path output_path;
  fs::path temporary;
  try {
    for (int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if (++index >= argc) throw std::runtime_error("Missing option value");
      const fs::path value = argv[index];
      if (option == "--weighted-ledger") {
        input_path = value;
      } else if (option == "--output") {
        output_path = value;
      } else {
        throw std::runtime_error("Unknown option: " + option);
      }
    }
    if (input_path.empty() || output_path.empty()) {
      throw std::runtime_error("Required: --weighted-ledger PATH --output PATH");
    }
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(input_path) ||
        fs::file_size(input_path) % sizeof(WeightedRecord)) {
      throw std::runtime_error("The weighted ledger has the wrong size");
    }
    if (fs::exists(output_path)) {
      throw std::runtime_error("Refusing to overwrite a projected ledger");
    }
    std::ifstream input(input_path, std::ios::binary);
    temporary = output_path.string() + ".tmp";
    if (fs::exists(temporary)) {
      throw std::runtime_error("A temporary projected ledger exists");
    }
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!input || !output) throw std::runtime_error("Could not open projection streams");
    WeightedRecord record{};
    WeightedRecord prior{};
    bool have_prior = false;
    std::uint64_t record_count = 0;
    std::uint64_t member_count = 0;
    while (input.read(reinterpret_cast<char*>(&record), sizeof(record))) {
      if (record.members == 0 || (have_prior && !less(prior, record)) ||
          std::numeric_limits<std::uint64_t>::max() - member_count <
              record.members) {
        throw std::runtime_error("A weighted projection input is invalid");
      }
      const ProjectedRecord projected{record.signature, record.mask};
      output.write(reinterpret_cast<const char*>(&projected), sizeof(projected));
      if (!output) throw std::runtime_error("Could not write a projected record");
      prior = record;
      have_prior = true;
      member_count += record.members;
      ++record_count;
    }
    if (!input.eof() || input.gcount() != 0) {
      throw std::runtime_error("The weighted ledger is truncated");
    }
    output.close();
    if (!output || fs::file_size(temporary) != record_count * sizeof(ProjectedRecord)) {
      throw std::runtime_error("The projected ledger has the wrong size");
    }
    fs::rename(temporary, output_path);
    temporary.clear();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"complete_length54_m09_weighted_ledger_projection_v1\",\n"
        << "  \"record_count\": " << record_count << ",\n"
        << "  \"source_member_count\": " << member_count << ",\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(input_path) << ",\n"
        << "  \"projected_ledger_bytes\": " << fs::file_size(output_path) << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::error_code ignored;
    if (!temporary.empty()) fs::remove(temporary, ignored);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
