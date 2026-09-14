#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <sys/resource.h>

namespace fs = std::filesystem;

namespace {

struct Mask {
  std::array<std::uint64_t, 8> words{};
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

struct RepeatedBucketRecord {
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t candidate_count = 0;
};

static_assert(sizeof(Record) == 96);
static_assert(sizeof(BucketRecord) == 176);
static_assert(sizeof(RepeatedBucketRecord) == 16);

struct Arguments {
  fs::path input_list;
  fs::path output;
  fs::path bucket_index;
  fs::path repeated_bucket_index;
  bool discard_bucket_index = false;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool record_less(const Record& left, const Record& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool record_equal(const Record& left, const Record& right) {
  return left.signature == right.signature && left.mask.words == right.mask.words;
}

struct HeapEntry {
  Record record{};
  std::size_t stream_index = 0;
};

struct HeapGreater {
  bool operator()(const HeapEntry& left, const HeapEntry& right) const {
    return record_less(right.record, left.record);
  }
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--discard-bucket-index") {
      result.discard_bucket_index = true;
      continue;
    }
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[++index];
    if (option == "--input-list") {
      result.input_list = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--bucket-index") {
      result.bucket_index = value;
    } else if (option == "--repeated-bucket-index") {
      result.repeated_bucket_index = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  const auto bucket_mode_count =
      static_cast<int>(!result.bucket_index.empty()) +
      static_cast<int>(!result.repeated_bucket_index.empty()) +
      static_cast<int>(result.discard_bucket_index);
  if (result.input_list.empty() || result.output.empty() ||
      bucket_mode_count != 1) {
    throw std::runtime_error(
        "Required: --input-list PATH --output PATH and exactly one of "
        "--bucket-index PATH, --repeated-bucket-index PATH, or "
        "--discard-bucket-index");
  }
  return result;
}

std::vector<fs::path> read_paths(const fs::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Could not open input-list file");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.empty()) throw std::runtime_error("The input list is empty");
  return result;
}

bool read_record(std::ifstream& stream, Record& record) {
  stream.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (stream.gcount() == 0 && stream.eof()) return false;
  if (stream.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("A signature shard is truncated");
  }
  return true;
}

struct FileDescriptorLimits {
  std::uint64_t soft_before = 0;
  std::uint64_t soft_after = 0;
};

std::uint64_t limit_value(rlim_t value) {
  return value == RLIM_INFINITY
      ? std::numeric_limits<std::uint64_t>::max()
      : static_cast<std::uint64_t>(value);
}

FileDescriptorLimits ensure_file_descriptor_capacity(
    std::size_t nonempty_shard_count) {
  constexpr rlim_t kDescriptorHeadroom = 32;
  struct rlimit limits {};
  if (getrlimit(RLIMIT_NOFILE, &limits) != 0) {
    throw std::runtime_error("Could not inspect the file-descriptor limit");
  }
  FileDescriptorLimits result;
  result.soft_before = limit_value(limits.rlim_cur);
  const auto required = static_cast<rlim_t>(nonempty_shard_count) +
                        kDescriptorHeadroom;
  if (limits.rlim_cur < required) {
    limits.rlim_cur = std::min(required, limits.rlim_max);
    if (setrlimit(RLIMIT_NOFILE, &limits) != 0) {
      throw std::runtime_error("Could not raise the file-descriptor limit");
    }
    if (getrlimit(RLIMIT_NOFILE, &limits) != 0) {
      throw std::runtime_error("Could not verify the file-descriptor limit");
    }
  }
  result.soft_after = limit_value(limits.rlim_cur);
  if (limits.rlim_cur < required) {
    throw std::runtime_error(
        "The hard file-descriptor limit is too low for a one-pass merge");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  fs::path output_temp;
  fs::path bucket_temp;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    const auto bucket_output = arguments.bucket_index.empty()
        ? arguments.repeated_bucket_index
        : arguments.bucket_index;
    if (fs::exists(arguments.output) ||
        (!arguments.discard_bucket_index && fs::exists(bucket_output))) {
      throw std::runtime_error("Refusing to overwrite a merged artifact");
    }
    const auto paths = read_paths(arguments.input_list);
    std::vector<std::uint64_t> shard_record_counts(paths.size());
    std::uint64_t input_record_count = 0;
    std::size_t nonempty_shard_count = 0;
    for (std::size_t index = 0; index < paths.size(); ++index) {
      if (!fs::is_regular_file(paths[index]) ||
          fs::file_size(paths[index]) % sizeof(Record) != 0) {
        throw std::runtime_error("A signature shard has the wrong byte length");
      }
      shard_record_counts[index] = fs::file_size(paths[index]) / sizeof(Record);
      input_record_count += shard_record_counts[index];
      if (shard_record_counts[index] != 0) ++nonempty_shard_count;
    }
    const auto descriptor_limits =
        ensure_file_descriptor_capacity(nonempty_shard_count);
    std::vector<std::ifstream> streams(paths.size());
    std::vector<Record> stream_priors(paths.size());
    std::vector<bool> stream_has_prior(paths.size(), false);
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapGreater> heap;
    for (std::size_t index = 0; index < paths.size(); ++index) {
      if (shard_record_counts[index] == 0) continue;
      streams[index].open(paths[index], std::ios::binary);
      if (!streams[index]) {
        throw std::runtime_error("Could not open a signature shard");
      }
      Record first;
      if (!read_record(streams[index], first)) {
        throw std::runtime_error("A nonempty signature shard became empty");
      }
      stream_priors[index] = first;
      stream_has_prior[index] = true;
      heap.push({first, index});
    }

    output_temp = arguments.output.string() + ".tmp";
    if (!arguments.discard_bucket_index) {
      bucket_temp = bucket_output.string() + ".tmp";
    }
    std::ofstream output(output_temp, std::ios::binary | std::ios::trunc);
    std::ofstream buckets;
    if (!arguments.discard_bucket_index) {
      buckets.open(bucket_temp, std::ios::binary | std::ios::trunc);
    }
    if (!output || (!arguments.discard_bucket_index && !buckets)) {
      throw std::runtime_error("Could not create merge outputs");
    }

    std::uint64_t unique_count = 0;
    std::uint64_t duplicate_count = 0;
    std::uint64_t consumed_count = 0;
    std::uint64_t bucket_count = 0;
    std::uint64_t repeated_bucket_count = 0;
    std::uint64_t largest_bucket = 0;
    bool has_global_prior = false;
    Record global_prior;
    bool has_bucket = false;
    BucketRecord bucket;

    const auto finish_bucket = [&]() {
      if (!has_bucket) return;
      const bool repeated = bucket.candidate_count > 1;
      if (repeated) ++repeated_bucket_count;
      if (!arguments.discard_bucket_index) {
        if (!arguments.bucket_index.empty()) {
          buckets.write(
              reinterpret_cast<const char*>(&bucket), sizeof(bucket));
        } else if (repeated) {
          const RepeatedBucketRecord repeated{
              bucket.ledger_offset_records, bucket.candidate_count};
          buckets.write(
              reinterpret_cast<const char*>(&repeated), sizeof(repeated));
        }
        if (!buckets) throw std::runtime_error("Could not write bucket index");
      }
      ++bucket_count;
      if (bucket.candidate_count > largest_bucket) {
        largest_bucket = bucket.candidate_count;
      }
      has_bucket = false;
    };

    while (!heap.empty()) {
      const auto entry = heap.top();
      heap.pop();
      ++consumed_count;
      if (has_global_prior && record_less(entry.record, global_prior)) {
        throw std::runtime_error("The global merge order regressed");
      }
      if (has_global_prior && record_equal(entry.record, global_prior)) {
        ++duplicate_count;
      } else {
        if (!has_bucket || bucket.signature != entry.record.signature) {
          finish_bucket();
          bucket = BucketRecord{};
          bucket.signature = entry.record.signature;
          bucket.ledger_offset_records = unique_count;
          bucket.candidate_count = 0;
          bucket.first_mask = entry.record.mask;
          has_bucket = true;
        }
        bucket.last_mask = entry.record.mask;
        ++bucket.candidate_count;
        output.write(reinterpret_cast<const char*>(&entry.record), sizeof(entry.record));
        if (!output) throw std::runtime_error("Could not write merged ledger");
        ++unique_count;
        global_prior = entry.record;
        has_global_prior = true;
      }

      Record next;
      auto& stream = streams[entry.stream_index];
      if (read_record(stream, next)) {
        if (stream_has_prior[entry.stream_index] &&
            !record_less(stream_priors[entry.stream_index], next)) {
          throw std::runtime_error("An input signature shard is not strictly sorted");
        }
        stream_priors[entry.stream_index] = next;
        stream_has_prior[entry.stream_index] = true;
        heap.push({next, entry.stream_index});
      }
    }
    finish_bucket();
    output.close();
    if (!arguments.discard_bucket_index) buckets.close();
    if (!output ||
        (!arguments.discard_bucket_index && !buckets) ||
        consumed_count != input_record_count ||
        fs::file_size(output_temp) != unique_count * sizeof(Record) ||
        (!arguments.discard_bucket_index &&
         fs::file_size(bucket_temp) !=
             (!arguments.bucket_index.empty()
                  ? bucket_count * sizeof(BucketRecord)
                  : repeated_bucket_count * sizeof(RepeatedBucketRecord))) ||
        unique_count + duplicate_count != input_record_count) {
      fs::remove(output_temp);
      if (!bucket_temp.empty()) fs::remove(bucket_temp);
      throw std::runtime_error("The merge accounting did not close");
    }
    fs::rename(output_temp, arguments.output);
    if (!arguments.discard_bucket_index) {
      fs::rename(bucket_temp, bucket_output);
    }
    const auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m09_signature_merge_v1\",\n"
              << "  \"input_shard_count\": " << paths.size() << ",\n"
              << "  \"nonempty_input_shard_count\": "
              << nonempty_shard_count << ",\n"
              << "  \"file_descriptor_soft_limit_before\": "
              << descriptor_limits.soft_before << ",\n"
              << "  \"file_descriptor_soft_limit_after\": "
              << descriptor_limits.soft_after << ",\n"
              << "  \"bucket_index_persisted\": "
              << (arguments.discard_bucket_index ? "false" : "true")
              << ",\n"
              << "  \"bucket_index_mode\": \""
              << (arguments.discard_bucket_index
                      ? "discarded"
                      : (!arguments.bucket_index.empty() ? "full"
                                                         : "repeated_only"))
              << "\",\n"
              << "  \"input_record_count\": " << input_record_count << ",\n"
              << "  \"unique_candidate_count\": " << unique_count << ",\n"
              << "  \"literal_duplicate_count\": " << duplicate_count << ",\n"
              << "  \"signature_bucket_count\": " << bucket_count << ",\n"
              << "  \"repeated_signature_bucket_count\": "
              << repeated_bucket_count << ",\n"
              << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::error_code ignored;
    if (!output_temp.empty()) fs::remove(output_temp, ignored);
    if (!bucket_temp.empty()) fs::remove(bucket_temp, ignored);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
