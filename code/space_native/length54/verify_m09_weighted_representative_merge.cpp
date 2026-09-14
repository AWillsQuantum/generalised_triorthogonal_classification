#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

using Signature = std::array<std::uint64_t, 4>;
using Mask = std::array<std::uint64_t, 8>;

struct WeightedRecord {
  Signature signature{};
  Mask mask{};
  std::uint64_t members = 0;
};

struct RepeatedRecord {
  std::uint64_t offset = 0;
  std::uint64_t count = 0;
};

static_assert(sizeof(WeightedRecord) == 104);
static_assert(sizeof(RepeatedRecord) == 16);

struct Arguments {
  fs::path inputs;
  fs::path merged;
  fs::path repeated;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[index];
    if (option == "--input-list") {
      result.inputs = value;
    } else if (option == "--merged-ledger") {
      result.merged = value;
    } else if (option == "--repeated-bucket-index") {
      result.repeated = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.inputs.empty() || result.merged.empty() || result.repeated.empty()) {
    throw std::runtime_error(
        "Required: --input-list PATH --merged-ledger PATH "
        "--repeated-bucket-index PATH");
  }
  return result;
}

std::vector<fs::path> paths(const fs::path& list) {
  std::ifstream input(list);
  if (!input) throw std::runtime_error("Could not open the input list");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.size() < 2) throw std::runtime_error("Too few weighted inputs");
  return result;
}

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

bool same_key(const WeightedRecord& left, const WeightedRecord& right) {
  return left.signature == right.signature && left.mask == right.mask;
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

struct Entry {
  WeightedRecord record{};
  std::size_t stream = 0;
};

struct Greater {
  bool operator()(const Entry& left, const Entry& right) const {
    return less(right.record, left.record);
  }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto input_paths = paths(arguments.inputs);
    if (!fs::is_regular_file(arguments.merged) ||
        fs::file_size(arguments.merged) % sizeof(WeightedRecord) ||
        !fs::is_regular_file(arguments.repeated) ||
        fs::file_size(arguments.repeated) % sizeof(RepeatedRecord)) {
      throw std::runtime_error("A weighted merge audit input has the wrong size");
    }
    std::vector<std::ifstream> inputs(input_paths.size());
    std::vector<WeightedRecord> priors(input_paths.size());
    std::vector<bool> have_prior(input_paths.size(), false);
    std::priority_queue<Entry, std::vector<Entry>, Greater> heap;
    for (std::size_t index = 0; index < input_paths.size(); ++index) {
      if (!fs::is_regular_file(input_paths[index]) ||
          fs::file_size(input_paths[index]) % sizeof(WeightedRecord)) {
        throw std::runtime_error("A weighted source has the wrong size");
      }
      inputs[index].open(input_paths[index], std::ios::binary);
      if (!inputs[index]) throw std::runtime_error("Could not open a source");
      WeightedRecord first;
      if (read_record(inputs[index], first, "A weighted source")) {
        if (first.members == 0) throw std::runtime_error("A source weight is zero");
        priors[index] = first;
        have_prior[index] = true;
        heap.push({first, index});
      }
    }
    std::ifstream merged(arguments.merged, std::ios::binary);
    std::ifstream repeated(arguments.repeated, std::ios::binary);
    if (!merged || !repeated) throw std::runtime_error("Could not open audit outputs");

    WeightedRecord pending{};
    bool have_pending = false;
    Signature bucket_signature{};
    bool have_bucket = false;
    std::uint64_t bucket_offset = 0;
    std::uint64_t bucket_size = 0;
    std::uint64_t input_records = 0;
    std::uint64_t input_members = 0;
    std::uint64_t unique_records = 0;
    std::uint64_t duplicate_records = 0;
    std::uint64_t output_members = 0;
    std::uint64_t buckets = 0;
    std::uint64_t repeated_buckets = 0;
    std::uint64_t repeated_representatives = 0;
    std::uint64_t singleton_buckets = 0;

    const auto close_bucket = [&]() {
      if (!have_bucket) return;
      ++buckets;
      if (bucket_size == 1) {
        ++singleton_buckets;
        return;
      }
      RepeatedRecord actual;
      if (!read_record(repeated, actual, "The repeated index") ||
          actual.offset != bucket_offset || actual.count != bucket_size) {
        throw std::runtime_error("A repeated-index record is incorrect");
      }
      ++repeated_buckets;
      repeated_representatives += bucket_size;
    };

    const auto check_pending = [&]() {
      if (!have_pending) return;
      WeightedRecord actual;
      if (!read_record(merged, actual, "The merged weighted ledger") ||
          actual.signature != pending.signature || actual.mask != pending.mask ||
          actual.members != pending.members) {
        throw std::runtime_error("A merged weighted record is incorrect");
      }
      if (!have_bucket || bucket_signature != pending.signature) {
        close_bucket();
        bucket_signature = pending.signature;
        bucket_offset = unique_records;
        bucket_size = 0;
        have_bucket = true;
      }
      ++unique_records;
      ++bucket_size;
      output_members += pending.members;
    };

    while (!heap.empty()) {
      auto entry = heap.top();
      heap.pop();
      ++input_records;
      if (entry.record.members == 0 ||
          std::numeric_limits<std::uint64_t>::max() - input_members <
              entry.record.members) {
        throw std::runtime_error("An input member count is invalid");
      }
      input_members += entry.record.members;
      if (have_pending && same_key(pending, entry.record)) {
        if (std::numeric_limits<std::uint64_t>::max() - pending.members <
            entry.record.members) {
          throw std::runtime_error("A merged member count overflows");
        }
        pending.members += entry.record.members;
        ++duplicate_records;
      } else {
        if (have_pending && less(entry.record, pending)) {
          throw std::runtime_error("The reconstructed merge order regressed");
        }
        check_pending();
        pending = entry.record;
        have_pending = true;
      }
      WeightedRecord next;
      if (read_record(inputs[entry.stream], next, "A weighted source")) {
        if (next.members == 0 ||
            (have_prior[entry.stream] && !less(priors[entry.stream], next))) {
          throw std::runtime_error("A weighted source is not strictly sorted");
        }
        priors[entry.stream] = next;
        have_prior[entry.stream] = true;
        heap.push({next, entry.stream});
      }
    }
    check_pending();
    close_bucket();
    WeightedRecord trailing_weighted;
    RepeatedRecord trailing_repeated;
    if (read_record(merged, trailing_weighted, "The merged weighted ledger") ||
        read_record(repeated, trailing_repeated, "The repeated index") ||
        input_records != unique_records + duplicate_records ||
        input_members != output_members ||
        unique_records != repeated_representatives + singleton_buckets ||
        fs::file_size(arguments.merged) != unique_records * sizeof(WeightedRecord) ||
        fs::file_size(arguments.repeated) !=
            repeated_buckets * sizeof(RepeatedRecord)) {
      throw std::runtime_error("The independent weighted merge audit did not close");
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_weighted_representative_merge_v1\",\n"
        << "  \"input_ledger_count\": " << input_paths.size() << ",\n"
        << "  \"input_record_count\": " << input_records << ",\n"
        << "  \"input_member_count\": " << input_members << ",\n"
        << "  \"unique_representative_count\": " << unique_records << ",\n"
        << "  \"literal_duplicate_count\": " << duplicate_records << ",\n"
        << "  \"output_member_count\": " << output_members << ",\n"
        << "  \"signature_bucket_count\": " << buckets << ",\n"
        << "  \"repeated_signature_bucket_count\": " << repeated_buckets << ",\n"
        << "  \"repeated_representative_count\": "
        << repeated_representatives << ",\n"
        << "  \"singleton_signature_bucket_count\": " << singleton_buckets
        << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
