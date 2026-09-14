#include <algorithm>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

struct RepeatedBucketRecord {
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t candidate_count = 0;
};

struct RangeRecord {
  std::uint64_t first_bucket = 0;
  std::uint64_t last_bucket = 0;
  std::uint64_t bucket_count = 0;
  std::uint64_t candidate_count = 0;
};

static_assert(sizeof(RepeatedBucketRecord) == 16);
static_assert(sizeof(RangeRecord) == 32);

struct Arguments {
  fs::path repeated_bucket_index;
  std::uint64_t ledger_record_count = 0;
  std::uint64_t target_candidates = 0;
  std::uint64_t maximum_buckets = 0;
  fs::path output;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--repeated-bucket-index") {
      result.repeated_bucket_index = value;
    } else if (option == "--ledger-record-count") {
      result.ledger_record_count = std::stoull(value);
    } else if (option == "--target-candidates") {
      result.target_candidates = std::stoull(value);
    } else if (option == "--maximum-buckets") {
      result.maximum_buckets = std::stoull(value);
    } else if (option == "--output") {
      result.output = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.repeated_bucket_index.empty() || result.ledger_record_count == 0 ||
      result.target_candidates == 0 || result.maximum_buckets == 0 ||
      result.output.empty()) {
    throw std::runtime_error(
        "Required: --repeated-bucket-index PATH --ledger-record-count N "
        "--target-candidates N --maximum-buckets N --output PATH");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  fs::path temporary;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.repeated_bucket_index) ||
        fs::file_size(arguments.repeated_bucket_index) %
                sizeof(RepeatedBucketRecord) !=
            0) {
      throw std::runtime_error("The repeated-bucket index has the wrong length");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite a quotient range plan");
    }
    temporary = arguments.output.string() + ".tmp";
    if (fs::exists(temporary)) {
      throw std::runtime_error("A temporary quotient range plan already exists");
    }
    std::ifstream input(arguments.repeated_bucket_index, std::ios::binary);
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!input || !output) throw std::runtime_error("Could not open a planner stream");

    const auto expected_bucket_count =
        fs::file_size(arguments.repeated_bucket_index) /
        sizeof(RepeatedBucketRecord);
    std::uint64_t bucket_count = 0;
    std::uint64_t candidate_count = 0;
    std::uint64_t range_count = 0;
    std::uint64_t range_first = 0;
    std::uint64_t range_buckets = 0;
    std::uint64_t range_candidates = 0;
    std::uint64_t previous_end = 0;
    std::uint64_t largest_bucket = 0;
    std::uint64_t largest_range_buckets = 0;
    std::uint64_t largest_range_candidates = 0;

    const auto flush = [&](std::uint64_t last_bucket) {
      if (range_buckets == 0) return;
      const RangeRecord range{
          range_first, last_bucket, range_buckets, range_candidates};
      output.write(reinterpret_cast<const char*>(&range), sizeof(range));
      if (!output) throw std::runtime_error("Could not write a quotient range");
      ++range_count;
      largest_range_buckets = std::max(largest_range_buckets, range_buckets);
      largest_range_candidates =
          std::max(largest_range_candidates, range_candidates);
      range_first = last_bucket;
      range_buckets = 0;
      range_candidates = 0;
    };

    RepeatedBucketRecord bucket;
    while (input.read(reinterpret_cast<char*>(&bucket), sizeof(bucket))) {
      if (bucket.candidate_count <= 1 ||
          bucket.ledger_offset_records < previous_end ||
          bucket.ledger_offset_records + bucket.candidate_count >
              arguments.ledger_record_count) {
        throw std::runtime_error("The repeated-bucket index is invalid");
      }
      previous_end = bucket.ledger_offset_records + bucket.candidate_count;
      if (range_buckets != 0 &&
          (range_buckets >= arguments.maximum_buckets ||
           range_candidates + bucket.candidate_count >
               arguments.target_candidates)) {
        flush(bucket_count);
      }
      ++bucket_count;
      ++range_buckets;
      candidate_count += bucket.candidate_count;
      range_candidates += bucket.candidate_count;
      largest_bucket = std::max(largest_bucket, bucket.candidate_count);
    }
    if (input.gcount() != 0 || !input.eof()) {
      throw std::runtime_error("The repeated-bucket index is truncated");
    }
    flush(bucket_count);
    output.close();
    if (!output || bucket_count != expected_bucket_count ||
        fs::file_size(temporary) != range_count * sizeof(RangeRecord)) {
      throw std::runtime_error("The quotient range-plan accounting did not close");
    }
    fs::rename(temporary, arguments.output);
    temporary.clear();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"complete_length54_m09_repeated_quotient_range_plan_v1\",\n"
        << "  \"repeated_bucket_count\": " << bucket_count << ",\n"
        << "  \"repeated_candidate_count\": " << candidate_count << ",\n"
        << "  \"range_count\": " << range_count << ",\n"
        << "  \"largest_bucket\": " << largest_bucket << ",\n"
        << "  \"largest_range_bucket_count\": " << largest_range_buckets
        << ",\n"
        << "  \"largest_range_candidate_count\": "
        << largest_range_candidates << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << ",\n"
        << "  \"buckets_per_second\": "
        << (elapsed == 0 ? 0 : bucket_count / elapsed) << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::error_code ignored;
    if (!temporary.empty()) fs::remove(temporary, ignored);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
