#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

struct Mask {
  std::array<std::uint64_t, 256> words{};
};

struct Record {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

struct BucketRecord {
  std::array<std::uint64_t, 4> signature{};
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t candidate_count = 0;
  Mask first_mask{};
  Mask last_mask{};
};

static_assert(sizeof(Record) == 2080);
static_assert(sizeof(BucketRecord) == 4144);

struct Arguments {
  fs::path ledger;
  fs::path bucket_index;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 255; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
  return result;
}

template <typename T>
bool read_record(std::ifstream& input, T& record) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("A merged signature artifact is truncated");
  }
  return true;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string option = argv[index];
    const fs::path value = argv[++index];
    if (option == "--ledger") {
      result.ledger = value;
    } else if (option == "--bucket-index") {
      result.bucket_index = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.ledger.empty() || result.bucket_index.empty()) {
    throw std::runtime_error("Required: --ledger PATH --bucket-index PATH");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.ledger) ||
        !fs::is_regular_file(arguments.bucket_index) ||
        fs::file_size(arguments.ledger) % sizeof(Record) != 0 ||
        fs::file_size(arguments.bucket_index) % sizeof(BucketRecord) != 0) {
      throw std::runtime_error("A merged signature artifact has a bad size");
    }
    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream buckets(arguments.bucket_index, std::ios::binary);
    if (!ledger || !buckets) throw std::runtime_error("Could not open merge artifacts");

    std::uint64_t candidate_count = 0;
    std::uint64_t bucket_count = 0;
    std::uint64_t largest_bucket = 0;
    bool has_prior_signature = false;
    std::array<std::uint64_t, 4> prior_signature{};
    BucketRecord bucket;
    while (read_record(buckets, bucket)) {
      if (bucket.ledger_offset_records != candidate_count ||
          bucket.candidate_count == 0 ||
          (has_prior_signature && !(prior_signature < bucket.signature)) ||
          mask_less(bucket.last_mask, bucket.first_mask)) {
        throw std::runtime_error("The bucket index order or accounting failed");
      }
      bool has_prior_mask = false;
      Mask prior_mask;
      for (std::uint64_t index = 0; index < bucket.candidate_count; ++index) {
        Record record;
        if (!read_record(ledger, record)) {
          throw std::runtime_error("The bucket index extends beyond the ledger");
        }
        if (record.signature != bucket.signature || mask_weight(record.mask) != 54 ||
            (has_prior_mask && !mask_less(prior_mask, record.mask)) ||
            (index == 0 && record.mask.words != bucket.first_mask.words) ||
            (index + 1 == bucket.candidate_count &&
             record.mask.words != bucket.last_mask.words)) {
          throw std::runtime_error("A ledger record does not match its bucket");
        }
        prior_mask = record.mask;
        has_prior_mask = true;
      }
      candidate_count += bucket.candidate_count;
      largest_bucket = std::max(largest_bucket, bucket.candidate_count);
      prior_signature = bucket.signature;
      has_prior_signature = true;
      ++bucket_count;
    }
    Record extra;
    if (read_record(ledger, extra)) {
      throw std::runtime_error("The ledger has records beyond the bucket index");
    }
    if (candidate_count * sizeof(Record) != fs::file_size(arguments.ledger) ||
        bucket_count * sizeof(BucketRecord) !=
            fs::file_size(arguments.bucket_index)) {
      throw std::runtime_error("The final merged sizes do not close");
    }
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    std::cout << "{\n"
              << "  \"status\": \"verified_length54_m14_merged_signature_ledger_v1\",\n"
              << "  \"candidate_count\": " << candidate_count << ",\n"
              << "  \"signature_bucket_count\": " << bucket_count << ",\n"
              << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "verify_m14_merged_signature_ledger: " << error.what() << '\n';
    return 1;
  }
}
