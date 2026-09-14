#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 8;
constexpr std::size_t kLedgerRecordBytes = 64;
constexpr std::size_t kAssignmentRecordBytes = 16;

using Signature = std::array<std::uint64_t, 4>;
using Mask = std::array<std::uint64_t, 4>;

struct LedgerRecord {
  Signature signature{};
  Mask mask{};
};

struct RepresentativeRecord {
  Signature signature{};
  Mask mask{};
  std::uint64_t member_count = 0;
};

static_assert(sizeof(LedgerRecord) == kLedgerRecordBytes);
static_assert(sizeof(RepresentativeRecord) == 72);

struct Arguments {
  fs::path ledger;
  fs::path assignments;
  fs::path representatives;
  std::uint64_t first_record = 0;
  std::uint64_t record_count = 0;
  std::uint64_t expected_buckets = 0;
  std::uint64_t expected_classes = 0;
  int support_size = 0;
  Signature expected_first_signature{};
  Signature expected_last_signature{};
  bool have_first_signature = false;
  bool have_last_signature = false;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 3; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool record_less(const LedgerRecord& left, const LedgerRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

void mask_insert(Mask& mask, std::uint8_t point) {
  mask[point >> 6] |= 1ULL << (point & 63);
}

bool mask_contains(const Mask& mask, std::uint8_t point) {
  return ((mask[point >> 6] >> (point & 63)) & 1ULL) != 0;
}

std::uint64_t mask_weight(const Mask& mask) {
  std::uint64_t result = 0;
  for (const auto word : mask) result += __builtin_popcountll(word);
  return result;
}

Signature parse_signature(const std::string& value) {
  if (value.size() != 64) {
    throw std::runtime_error("A signature must contain exactly 64 hex digits");
  }
  Signature result{};
  for (int index = 0; index < 4; ++index) {
    result[index] = std::stoull(value.substr(16 * index, 16), nullptr, 16);
  }
  return result;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--ledger") {
      result.ledger = value;
    } else if (option == "--assignments") {
      result.assignments = value;
    } else if (option == "--representatives") {
      result.representatives = value;
    } else if (option == "--first-record") {
      result.first_record = std::stoull(value);
    } else if (option == "--record-count") {
      result.record_count = std::stoull(value);
    } else if (option == "--expected-buckets") {
      result.expected_buckets = std::stoull(value);
    } else if (option == "--expected-classes") {
      result.expected_classes = std::stoull(value);
    } else if (option == "--support-size") {
      result.support_size = std::stoi(value);
    } else if (option == "--expected-first-signature") {
      result.expected_first_signature = parse_signature(value);
      result.have_first_signature = true;
    } else if (option == "--expected-last-signature") {
      result.expected_last_signature = parse_signature(value);
      result.have_last_signature = true;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.ledger.empty() || result.assignments.empty() ||
      result.representatives.empty() || result.record_count == 0 ||
      result.expected_buckets == 0 || result.expected_classes == 0 ||
      result.support_size <= 0 || result.support_size > 256 ||
      !result.have_first_signature || !result.have_last_signature) {
    throw std::runtime_error(
        "Required: --ledger PATH --assignments PATH --representatives PATH "
        "--first-record N --record-count N --expected-buckets N "
        "--expected-classes N --support-size N "
        "--expected-first-signature HEX --expected-last-signature HEX");
  }
  return result;
}

bool insert_independent(
    std::array<std::uint8_t, kDimension>& basis, std::uint8_t vector) {
  auto residual = vector;
  for (int bit = kDimension - 1; bit >= 0; --bit) {
    if (((residual >> bit) & 1U) == 0) continue;
    if (basis[bit] != 0) {
      residual ^= basis[bit];
    } else {
      basis[bit] = residual;
      return true;
    }
  }
  return false;
}

bool linear_images_are_invertible(const std::array<std::uint8_t, 8>& images) {
  std::array<std::uint8_t, kDimension> basis{};
  for (const auto image : images) {
    if (!insert_independent(basis, image)) return false;
  }
  return true;
}

Mask transform_mask(
    const Mask& source,
    std::uint8_t translation,
    const std::array<std::uint8_t, 8>& images) {
  Mask result{};
  for (int point = 0; point < 256; ++point) {
    if (!mask_contains(source, static_cast<std::uint8_t>(point))) continue;
    auto transformed = translation;
    for (int bit = 0; bit < kDimension; ++bit) {
      if (((point >> bit) & 1) != 0) transformed ^= images[bit];
    }
    if (mask_contains(result, transformed)) {
      throw std::runtime_error("An affine witness is not injective on its support");
    }
    mask_insert(result, transformed);
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.ledger) ||
        !fs::is_regular_file(arguments.assignments)) {
      throw std::runtime_error("A quotient verification input does not exist");
    }
    if (fs::exists(arguments.representatives)) {
      throw std::runtime_error("Refusing to overwrite representative output");
    }
    const auto ledger_records = fs::file_size(arguments.ledger) / kLedgerRecordBytes;
    if (fs::file_size(arguments.ledger) % kLedgerRecordBytes != 0 ||
        arguments.first_record + arguments.record_count > ledger_records ||
        fs::file_size(arguments.assignments) !=
            arguments.record_count * kAssignmentRecordBytes) {
      throw std::runtime_error("A quotient verification input has the wrong size");
    }

    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream assignments(arguments.assignments, std::ios::binary);
    ledger.seekg(static_cast<std::streamoff>(
        arguments.first_record * kLedgerRecordBytes));
    const auto temporary = arguments.representatives.string() + ".tmp";
    std::ofstream representatives(temporary, std::ios::binary | std::ios::trunc);
    if (!ledger || !assignments || !representatives) {
      throw std::runtime_error("Could not open quotient verification files");
    }

    Signature current_signature{};
    Signature first_signature{};
    Signature last_signature{};
    bool have_signature = false;
    LedgerRecord prior{};
    bool have_prior = false;
    std::vector<Mask> class_representatives;
    std::vector<std::uint64_t> class_member_counts;
    std::map<std::uint64_t, std::uint64_t> class_count_distribution;
    std::uint64_t bucket_count = 0;
    std::uint64_t class_count = 0;
    std::uint64_t member_count = 0;

    const auto close_bucket = [&]() {
      if (!have_signature) return;
      ++class_count_distribution[class_representatives.size()];
      for (std::size_t index = 0; index < class_representatives.size(); ++index) {
        RepresentativeRecord record{
            current_signature,
            class_representatives[index],
            class_member_counts[index]};
        representatives.write(
            reinterpret_cast<const char*>(&record), sizeof(record));
        if (!representatives) {
          throw std::runtime_error("Could not write a verified representative");
        }
        ++class_count;
        member_count += class_member_counts[index];
      }
    };

    for (std::uint64_t index = 0; index < arguments.record_count; ++index) {
      LedgerRecord ledger_record;
      std::array<std::uint8_t, kAssignmentRecordBytes> assignment{};
      ledger.read(reinterpret_cast<char*>(&ledger_record), sizeof(ledger_record));
      assignments.read(
          reinterpret_cast<char*>(assignment.data()), assignment.size());
      if (!ledger || !assignments) {
        throw std::runtime_error("A quotient verification input is truncated");
      }
      if (have_prior && !record_less(prior, ledger_record)) {
        throw std::runtime_error("The quotient ledger is not strictly sorted");
      }
      prior = ledger_record;
      have_prior = true;
      if (!have_signature || ledger_record.signature != current_signature) {
        close_bucket();
        current_signature = ledger_record.signature;
        class_representatives.clear();
        class_member_counts.clear();
        ++bucket_count;
        if (!have_signature) first_signature = ledger_record.signature;
        last_signature = ledger_record.signature;
        have_signature = true;
      }

      const std::uint16_t class_index = static_cast<std::uint16_t>(
          assignment[0] | (static_cast<std::uint16_t>(assignment[1]) << 8));
      if (class_index > class_representatives.size()) {
        throw std::runtime_error("A quotient assignment skips a class index");
      }
      if (class_index == class_representatives.size()) {
        class_representatives.push_back(ledger_record.mask);
        class_member_counts.push_back(0);
      }
      std::array<std::uint8_t, 8> images{};
      for (int bit = 0; bit < 8; ++bit) images[bit] = assignment[3 + bit];
      if (!linear_images_are_invertible(images) ||
          mask_weight(ledger_record.mask) !=
              static_cast<std::uint64_t>(arguments.support_size) ||
          transform_mask(ledger_record.mask, assignment[2], images) !=
              class_representatives[class_index] ||
          std::any_of(
              assignment.begin() + 11,
              assignment.end(),
              [](std::uint8_t value) { return value != 0; })) {
        throw std::runtime_error("A quotient assignment witness is invalid");
      }
      ++class_member_counts[class_index];
    }
    close_bucket();
    representatives.close();
    if (!representatives || bucket_count != arguments.expected_buckets ||
        class_count != arguments.expected_classes ||
        member_count != arguments.record_count ||
        first_signature != arguments.expected_first_signature ||
        last_signature != arguments.expected_last_signature ||
        fs::file_size(temporary) != class_count * sizeof(RepresentativeRecord)) {
      fs::remove(temporary);
      throw std::runtime_error("The quotient verification accounting did not close");
    }
    fs::rename(temporary, arguments.representatives);
    const auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::cout << "{\n"
              << "  \"status\": \"verified_m08_exact_quotient_batch_v1\",\n"
              << "  \"candidate_count\": " << arguments.record_count << ",\n"
              << "  \"bucket_count\": " << bucket_count << ",\n"
              << "  \"affine_class_count\": " << class_count << ",\n"
              << "  \"verified_witness_count\": " << arguments.record_count
              << ",\n"
              << "  \"representative_member_count\": " << member_count << ",\n"
              << "  \"bucket_affine_class_count_distribution\": {";
    bool first = true;
    for (const auto& [classes, frequency] : class_count_distribution) {
      if (!first) std::cout << ',';
      std::cout << "\n    \"" << classes << "\": " << frequency;
      first = false;
    }
    if (!first) std::cout << '\n';
    std::cout << "  },\n"
              << "  \"representative_record_bytes\": 72,\n"
              << "  \"representative_bytes\": "
              << fs::file_size(arguments.representatives) << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
