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
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Mask {
  std::array<std::uint64_t, 8> words{};
};

struct WeightedRecord {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
  std::uint64_t member_count = 0;
};

struct RepeatedBucketRecord {
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t representative_count = 0;
};

static_assert(sizeof(WeightedRecord) == 104);
static_assert(sizeof(RepeatedBucketRecord) == 16);

struct Arguments {
  fs::path input_list;
  fs::path output;
  fs::path repeated_bucket_index;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[++index];
    if (option == "--input-list") {
      result.input_list = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--repeated-bucket-index") {
      result.repeated_bucket_index = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.input_list.empty() || result.output.empty() ||
      result.repeated_bucket_index.empty()) {
    throw std::runtime_error(
        "Required: --input-list PATH --output PATH "
        "--repeated-bucket-index PATH");
  }
  return result;
}

std::vector<fs::path> read_paths(const fs::path& list) {
  std::ifstream input(list);
  if (!input) throw std::runtime_error("Could not open the weighted input list");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.size() < 2) {
    throw std::runtime_error("At least two weighted ledgers are required");
  }
  return result;
}

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool record_less(const WeightedRecord& left, const WeightedRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool same_key(const WeightedRecord& left, const WeightedRecord& right) {
  return left.signature == right.signature && left.mask.words == right.mask.words;
}

bool read_record(std::ifstream& input, WeightedRecord& record) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("A weighted input ledger is truncated");
  }
  if (record.member_count == 0) {
    throw std::runtime_error("A weighted input record has zero members");
  }
  return true;
}

struct HeapEntry {
  WeightedRecord record{};
  std::size_t stream_index = 0;
};

struct HeapGreater {
  bool operator()(const HeapEntry& left, const HeapEntry& right) const {
    return record_less(right.record, left.record);
  }
};

}  // namespace

int main(int argc, char** argv) {
  fs::path output_temporary;
  fs::path index_temporary;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto paths = read_paths(arguments.input_list);
    if (fs::exists(arguments.output) ||
        fs::exists(arguments.repeated_bucket_index)) {
      throw std::runtime_error("Refusing to overwrite a weighted merge output");
    }
    std::vector<std::ifstream> streams(paths.size());
    std::vector<WeightedRecord> stream_priors(paths.size());
    std::vector<bool> stream_has_prior(paths.size(), false);
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapGreater> heap;
    std::uint64_t input_record_count = 0;
    for (std::size_t index = 0; index < paths.size(); ++index) {
      if (!fs::is_regular_file(paths[index]) ||
          fs::file_size(paths[index]) % sizeof(WeightedRecord)) {
        throw std::runtime_error("A weighted input has the wrong byte length");
      }
      streams[index].open(paths[index], std::ios::binary);
      if (!streams[index]) throw std::runtime_error("Could not open an input ledger");
      WeightedRecord first;
      if (read_record(streams[index], first)) {
        stream_priors[index] = first;
        stream_has_prior[index] = true;
        heap.push({first, index});
      }
      input_record_count += fs::file_size(paths[index]) / sizeof(WeightedRecord);
    }

    output_temporary = arguments.output.string() + ".tmp";
    index_temporary = arguments.repeated_bucket_index.string() + ".tmp";
    if (fs::exists(output_temporary) || fs::exists(index_temporary)) {
      throw std::runtime_error("A temporary weighted merge output exists");
    }
    std::ofstream output(output_temporary, std::ios::binary | std::ios::trunc);
    std::ofstream repeated(index_temporary, std::ios::binary | std::ios::trunc);
    if (!output || !repeated) throw std::runtime_error("Could not open merge outputs");

    WeightedRecord pending{};
    bool have_pending = false;
    std::array<std::uint64_t, 4> bucket_signature{};
    bool have_bucket = false;
    std::uint64_t bucket_offset = 0;
    std::uint64_t bucket_size = 0;
    std::uint64_t consumed_count = 0;
    std::uint64_t input_member_count = 0;
    std::uint64_t unique_count = 0;
    std::uint64_t output_member_count = 0;
    std::uint64_t literal_duplicate_count = 0;
    std::uint64_t signature_bucket_count = 0;
    std::uint64_t repeated_bucket_count = 0;
    std::uint64_t repeated_representative_count = 0;
    std::uint64_t singleton_bucket_count = 0;
    std::uint64_t largest_bucket = 0;

    const auto finish_bucket = [&]() {
      if (!have_bucket) return;
      ++signature_bucket_count;
      largest_bucket = std::max(largest_bucket, bucket_size);
      if (bucket_size > 1) {
        const RepeatedBucketRecord record{bucket_offset, bucket_size};
        repeated.write(reinterpret_cast<const char*>(&record), sizeof(record));
        if (!repeated) throw std::runtime_error("Could not write a repeated index");
        ++repeated_bucket_count;
        repeated_representative_count += bucket_size;
      } else {
        ++singleton_bucket_count;
      }
    };

    const auto finish_record = [&]() {
      if (!have_pending) return;
      if (!have_bucket || bucket_signature != pending.signature) {
        finish_bucket();
        bucket_signature = pending.signature;
        bucket_offset = unique_count;
        bucket_size = 0;
        have_bucket = true;
      }
      output.write(reinterpret_cast<const char*>(&pending), sizeof(pending));
      if (!output) throw std::runtime_error("Could not write the weighted merge");
      ++unique_count;
      ++bucket_size;
      output_member_count += pending.member_count;
    };

    while (!heap.empty()) {
      auto entry = heap.top();
      heap.pop();
      ++consumed_count;
      if (std::numeric_limits<std::uint64_t>::max() - input_member_count <
          entry.record.member_count) {
        throw std::runtime_error("Input member count overflow");
      }
      input_member_count += entry.record.member_count;
      if (have_pending && same_key(pending, entry.record)) {
        if (std::numeric_limits<std::uint64_t>::max() - pending.member_count <
            entry.record.member_count) {
          throw std::runtime_error("Merged member count overflow");
        }
        pending.member_count += entry.record.member_count;
        ++literal_duplicate_count;
      } else {
        if (have_pending && record_less(entry.record, pending)) {
          throw std::runtime_error("The weighted global merge order regressed");
        }
        finish_record();
        pending = entry.record;
        have_pending = true;
      }

      auto& stream = streams[entry.stream_index];
      WeightedRecord next;
      if (read_record(stream, next)) {
        if (stream_has_prior[entry.stream_index] &&
            !record_less(stream_priors[entry.stream_index], next)) {
          throw std::runtime_error("A weighted input is not strictly sorted");
        }
        stream_priors[entry.stream_index] = next;
        stream_has_prior[entry.stream_index] = true;
        heap.push({next, entry.stream_index});
      }
    }
    finish_record();
    finish_bucket();
    output.close();
    repeated.close();
    if (!output || !repeated || consumed_count != input_record_count ||
        input_record_count != unique_count + literal_duplicate_count ||
        input_member_count != output_member_count ||
        unique_count != repeated_representative_count + singleton_bucket_count ||
        fs::file_size(output_temporary) != unique_count * sizeof(WeightedRecord) ||
        fs::file_size(index_temporary) !=
            repeated_bucket_count * sizeof(RepeatedBucketRecord)) {
      throw std::runtime_error("The weighted merge accounting did not close");
    }
    fs::rename(output_temporary, arguments.output);
    output_temporary.clear();
    fs::rename(index_temporary, arguments.repeated_bucket_index);
    index_temporary.clear();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"complete_length54_m09_weighted_representative_merge_v1\",\n"
        << "  \"input_ledger_count\": " << paths.size() << ",\n"
        << "  \"input_record_count\": " << input_record_count << ",\n"
        << "  \"input_member_count\": " << input_member_count << ",\n"
        << "  \"unique_representative_count\": " << unique_count << ",\n"
        << "  \"literal_duplicate_count\": " << literal_duplicate_count << ",\n"
        << "  \"output_member_count\": " << output_member_count << ",\n"
        << "  \"signature_bucket_count\": " << signature_bucket_count << ",\n"
        << "  \"repeated_signature_bucket_count\": " << repeated_bucket_count << ",\n"
        << "  \"repeated_representative_count\": "
        << repeated_representative_count << ",\n"
        << "  \"singleton_signature_bucket_count\": "
        << singleton_bucket_count << ",\n"
        << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(arguments.output)
        << ",\n"
        << "  \"repeated_bucket_index_bytes\": "
        << fs::file_size(arguments.repeated_bucket_index) << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::error_code ignored;
    if (!output_temporary.empty()) fs::remove(output_temporary, ignored);
    if (!index_temporary.empty()) fs::remove(index_temporary, ignored);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
