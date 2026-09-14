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

namespace fs = std::filesystem;

namespace {

constexpr int kSupportSize = 54;

struct Mask {
  std::array<std::uint64_t, 8> words{};
};

struct Record {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

static_assert(sizeof(Mask) == 64);
static_assert(sizeof(Record) == 96);

struct Arguments {
  fs::path masks;
  fs::path signatures;
};

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

bool record_less(const Record& left, const Record& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
  return result;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string option = argv[index];
    const fs::path value = argv[++index];
    if (option == "--masks") {
      result.masks = value;
    } else if (option == "--signatures") {
      result.signatures = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.masks.empty() || result.signatures.empty()) {
    throw std::runtime_error("Required: --masks PATH --signatures PATH");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.masks) ||
        !fs::is_regular_file(arguments.signatures) ||
        fs::file_size(arguments.masks) % sizeof(Mask) != 0 ||
        fs::file_size(arguments.signatures) % sizeof(Record) != 0) {
      throw std::runtime_error("An input has an invalid byte count");
    }
    const std::uint64_t mask_count =
        fs::file_size(arguments.masks) / sizeof(Mask);
    const std::uint64_t signature_count =
        fs::file_size(arguments.signatures) / sizeof(Record);
    if (mask_count != signature_count) {
      throw std::runtime_error("Mask and signature record counts differ");
    }

    std::ifstream masks(arguments.masks, std::ios::binary);
    Mask prior_mask;
    bool has_prior_mask = false;
    for (std::uint64_t index = 0; index < mask_count; ++index) {
      Mask mask;
      masks.read(reinterpret_cast<char*>(&mask), sizeof(mask));
      if (!masks || mask_weight(mask) != kSupportSize ||
          (has_prior_mask && !mask_less(prior_mask, mask))) {
        throw std::runtime_error("The mask ledger is malformed or not strictly sorted");
      }
      prior_mask = mask;
      has_prior_mask = true;
    }
    if (masks.peek() != std::ifstream::traits_type::eof()) {
      throw std::runtime_error("The mask ledger has trailing data");
    }
    masks.close();

    std::vector<Mask> discovered;
    discovered.reserve(static_cast<std::size_t>(signature_count));
    std::ifstream signatures(arguments.signatures, std::ios::binary);
    Record prior_record;
    bool has_prior_record = false;
    for (std::uint64_t index = 0; index < signature_count; ++index) {
      Record record;
      signatures.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (!signatures || mask_weight(record.mask) != kSupportSize ||
          (has_prior_record && !record_less(prior_record, record))) {
        throw std::runtime_error(
            "The signature ledger is malformed or not strictly sorted");
      }
      discovered.push_back(record.mask);
      prior_record = record;
      has_prior_record = true;
    }
    if (signatures.peek() != std::ifstream::traits_type::eof()) {
      throw std::runtime_error("The signature ledger has trailing data");
    }
    signatures.close();

    std::sort(discovered.begin(), discovered.end(), mask_less);
    if (std::adjacent_find(
            discovered.begin(), discovered.end(), mask_equal) !=
        discovered.end()) {
      throw std::runtime_error("The signature ledger contains duplicate masks");
    }
    masks.open(arguments.masks, std::ios::binary);
    for (const auto& discovered_mask : discovered) {
      Mask expected;
      masks.read(reinterpret_cast<char*>(&expected), sizeof(expected));
      if (!masks || !mask_equal(expected, discovered_mask)) {
        throw std::runtime_error(
            "Signature masks do not exactly partition the input mask ledger");
      }
    }
    if (masks.peek() != std::ifstream::traits_type::eof()) {
      throw std::runtime_error("The comparison mask ledger has trailing data");
    }

    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_signature_batch_v1\",\n"
        << "  \"candidate_count\": " << mask_count << ",\n"
        << "  \"every_input_mask_matched_once\": true,\n"
        << "  \"mask_input_strictly_sorted\": true,\n"
        << "  \"signature_records_strictly_sorted\": true,\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "m09_verify_signature_batch: " << error.what() << "\n";
    return 1;
  }
}
