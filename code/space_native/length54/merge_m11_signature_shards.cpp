#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Mask {
  std::array<std::uint64_t, 32> words{};
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

static_assert(sizeof(Record) == 288);
static_assert(sizeof(BucketRecord) == 560);

struct Arguments {
  fs::path input_list;
  fs::path output;
  fs::path bucket_index;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 31; index >= 0; --index) {
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
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[++index];
    if (option == "--input-list") {
      result.input_list = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--bucket-index") {
      result.bucket_index = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.input_list.empty() || result.output.empty() ||
      result.bucket_index.empty()) {
    throw std::runtime_error(
        "Required: --input-list PATH --output PATH --bucket-index PATH");
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

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (fs::exists(arguments.output) || fs::exists(arguments.bucket_index)) {
      throw std::runtime_error("Refusing to overwrite a merged artifact");
    }
    const auto paths = read_paths(arguments.input_list);
    std::vector<std::ifstream> streams;
    streams.reserve(paths.size());
    std::vector<Record> stream_priors(paths.size());
    std::vector<bool> stream_has_prior(paths.size(), false);
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapGreater> heap;
    std::uint64_t input_record_count = 0;
    for (std::size_t index = 0; index < paths.size(); ++index) {
      if (!fs::is_regular_file(paths[index]) ||
          fs::file_size(paths[index]) % sizeof(Record) != 0) {
        throw std::runtime_error("A signature shard has the wrong byte length");
      }
      input_record_count += fs::file_size(paths[index]) / sizeof(Record);
      streams.emplace_back(paths[index], std::ios::binary);
      if (!streams.back()) throw std::runtime_error("Could not open a signature shard");
      Record first;
      if (read_record(streams.back(), first)) {
        stream_priors[index] = first;
        stream_has_prior[index] = true;
        heap.push({first, index});
      }
    }

    const auto output_temp = arguments.output.string() + ".tmp";
    const auto bucket_temp = arguments.bucket_index.string() + ".tmp";
    std::ofstream output(output_temp, std::ios::binary | std::ios::trunc);
    std::ofstream buckets(bucket_temp, std::ios::binary | std::ios::trunc);
    if (!output || !buckets) throw std::runtime_error("Could not create merge outputs");

    std::uint64_t unique_count = 0;
    std::uint64_t duplicate_count = 0;
    std::uint64_t consumed_count = 0;
    std::uint64_t bucket_count = 0;
    std::uint64_t largest_bucket = 0;
    bool has_global_prior = false;
    Record global_prior;
    bool has_bucket = false;
    BucketRecord bucket;

    const auto finish_bucket = [&]() {
      if (!has_bucket) return;
      buckets.write(reinterpret_cast<const char*>(&bucket), sizeof(bucket));
      if (!buckets) throw std::runtime_error("Could not write bucket index");
      ++bucket_count;
      largest_bucket = std::max(largest_bucket, bucket.candidate_count);
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
        heap.push({next, entry.stream_index});
      }
    }
    finish_bucket();
    output.close();
    buckets.close();
    if (!output || !buckets || consumed_count != input_record_count ||
        fs::file_size(output_temp) != unique_count * sizeof(Record) ||
        fs::file_size(bucket_temp) != bucket_count * sizeof(BucketRecord) ||
        unique_count + duplicate_count != input_record_count) {
      fs::remove(output_temp);
      fs::remove(bucket_temp);
      throw std::runtime_error("The merge accounting did not close");
    }
    fs::rename(output_temp, arguments.output);
    fs::rename(bucket_temp, arguments.bucket_index);
    const auto elapsed = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - started)
                             .count();
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m11_signature_merge_v1\",\n"
              << "  \"input_shard_count\": " << paths.size() << ",\n"
              << "  \"input_record_count\": " << input_record_count << ",\n"
              << "  \"unique_candidate_count\": " << unique_count << ",\n"
              << "  \"literal_duplicate_count\": " << duplicate_count << ",\n"
              << "  \"signature_bucket_count\": " << bucket_count << ",\n"
              << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
