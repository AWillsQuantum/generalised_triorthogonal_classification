#include <algorithm>
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

struct AssignmentRecord {
  std::uint32_t local_class = 0;
  std::uint16_t translation = 0;
  std::array<std::uint16_t, 9> basis_images{};
};

struct ClassRecord {
  std::uint64_t bucket = 0;
  std::uint32_t local_class = 0;
  std::uint32_t representative_index = 0;
  Mask representative{};
  std::uint64_t record_members = 0;
};

static_assert(sizeof(WeightedRecord) == 104);
static_assert(sizeof(RepeatedRecord) == 16);
static_assert(sizeof(AssignmentRecord) == 24);
static_assert(sizeof(ClassRecord) == 88);

struct Arguments {
  fs::path weighted_source;
  fs::path repeated;
  fs::path assignments;
  fs::path classes;
  fs::path output;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[index];
    if (option == "--weighted-source") result.weighted_source = value;
    else if (option == "--repeated-bucket-index") result.repeated = value;
    else if (option == "--assignment-input-list") result.assignments = value;
    else if (option == "--class-input-list") result.classes = value;
    else if (option == "--output") result.output = value;
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.weighted_source.empty() || result.repeated.empty() ||
      result.assignments.empty() || result.classes.empty() ||
      result.output.empty()) {
    throw std::runtime_error(
        "Required: --weighted-source PATH --repeated-bucket-index PATH "
        "--assignment-input-list PATH --class-input-list PATH --output PATH");
  }
  return result;
}

template <typename Record>
bool read_fixed(std::ifstream& input, Record& record, const char* label) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error(std::string(label) + " is truncated");
  }
  return true;
}

std::vector<fs::path> read_paths(const fs::path& list, std::size_t record_bytes,
                                 const char* label) {
  std::ifstream input(list);
  if (!input) throw std::runtime_error(std::string("Could not open ") + label);
  std::vector<fs::path> paths;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) paths.emplace_back(line);
  }
  if (paths.empty()) throw std::runtime_error(std::string(label) + " is empty");
  for (const auto& path : paths) {
    if (!fs::is_regular_file(path) || fs::file_size(path) % record_bytes) {
      throw std::runtime_error(std::string(label) + " has an invalid shard");
    }
  }
  return paths;
}

template <typename Record>
class RecordStream {
 public:
  explicit RecordStream(std::vector<fs::path> paths)
      : paths_(std::move(paths)) {}

  bool next(Record& record) {
    while (true) {
      if (!stream_.is_open()) {
        if (next_path_ == paths_.size()) return false;
        stream_.open(paths_[next_path_++], std::ios::binary);
        if (!stream_) throw std::runtime_error("Could not open an input shard");
      }
      stream_.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (stream_.gcount() == static_cast<std::streamsize>(sizeof(record))) {
        return true;
      }
      if (stream_.gcount() != 0 || !stream_.eof()) {
        throw std::runtime_error("An input shard is truncated");
      }
      stream_.close();
      stream_.clear();
    }
  }

 private:
  std::vector<fs::path> paths_;
  std::size_t next_path_ = 0;
  std::ifstream stream_;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool weighted_less(const WeightedRecord& left, const WeightedRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool identity(const AssignmentRecord& assignment) {
  if (assignment.translation != 0) return false;
  for (int bit = 0; bit < 9; ++bit) {
    if (assignment.basis_images[bit] != (1U << bit)) return false;
  }
  return true;
}

void checked_add(std::uint64_t& target, std::uint64_t value) {
  if (std::numeric_limits<std::uint64_t>::max() - target < value) {
    throw std::runtime_error("A member-count sum overflowed");
  }
  target += value;
}

void write_record(std::ofstream& output, const WeightedRecord& record) {
  output.write(reinterpret_cast<const char*>(&record), sizeof(record));
  if (!output) throw std::runtime_error("Could not write a representative");
}

}  // namespace

int main(int argc, char** argv) {
  fs::path temporary;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.weighted_source) ||
        fs::file_size(arguments.weighted_source) % sizeof(WeightedRecord) ||
        !fs::is_regular_file(arguments.repeated) ||
        fs::file_size(arguments.repeated) % sizeof(RepeatedRecord)) {
      throw std::runtime_error("A global representative input has the wrong size");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite a representative ledger");
    }
    const auto assignment_paths = read_paths(
        arguments.assignments, sizeof(AssignmentRecord), "the assignment list");
    const auto class_paths =
        read_paths(arguments.classes, sizeof(ClassRecord), "the class list");
    RecordStream<AssignmentRecord> assignments(assignment_paths);
    RecordStream<ClassRecord> classes(class_paths);
    std::ifstream source(arguments.weighted_source, std::ios::binary);
    std::ifstream repeated(arguments.repeated, std::ios::binary);
    temporary = arguments.output.string() + ".tmp";
    if (fs::exists(temporary)) {
      throw std::runtime_error("A temporary representative ledger exists");
    }
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!source || !repeated || !output) {
      throw std::runtime_error("Could not open a global representative stream");
    }

    WeightedRecord current{};
    bool have_current = read_fixed(source, current, "The weighted source");
    WeightedRecord prior{};
    bool have_prior = false;
    RepeatedRecord next_repeated{};
    bool have_repeated = read_fixed(repeated, next_repeated,
                                    "The repeated-bucket index");
    ClassRecord next_class{};
    bool have_class = classes.next(next_class);
    std::uint64_t source_records = 0;
    std::uint64_t source_members = 0;
    std::uint64_t signature_buckets = 0;
    std::uint64_t repeated_buckets = 0;
    std::uint64_t repeated_records = 0;
    std::uint64_t singleton_classes = 0;
    std::uint64_t exact_classes = 0;
    std::uint64_t output_members = 0;
    std::uint64_t split_buckets = 0;
    std::uint64_t assignment_records = 0;
    std::uint64_t largest_bucket = 0;

    while (have_current) {
      const auto bucket_offset = source_records;
      const auto signature = current.signature;
      std::vector<WeightedRecord> bucket;
      do {
        if (current.members == 0 ||
            (have_prior && !weighted_less(prior, current))) {
          throw std::runtime_error(
              "The weighted source is not strictly ordered and positive");
        }
        prior = current;
        have_prior = true;
        checked_add(source_members, current.members);
        bucket.push_back(current);
        ++source_records;
        have_current = read_fixed(source, current, "The weighted source");
      } while (have_current && current.signature == signature);
      ++signature_buckets;
      largest_bucket = std::max<std::uint64_t>(largest_bucket, bucket.size());
      if (have_repeated && next_repeated.offset < bucket_offset) {
        throw std::runtime_error("The repeated index is not ordered");
      }
      if (bucket.size() == 1) {
        if (have_repeated && next_repeated.offset == bucket_offset) {
          throw std::runtime_error("The repeated index contains a singleton");
        }
        write_record(output, bucket.front());
        checked_add(output_members, bucket.front().members);
        ++singleton_classes;
        continue;
      }
      if (!have_repeated || next_repeated.offset != bucket_offset ||
          next_repeated.count != bucket.size()) {
        throw std::runtime_error("The repeated index does not bind a source bucket");
      }

      std::vector<ClassRecord> bucket_classes;
      while (have_class && next_class.bucket == repeated_buckets) {
        if (next_class.local_class != bucket_classes.size() ||
            next_class.representative_index >= bucket.size() ||
            next_class.representative !=
                bucket[next_class.representative_index].mask ||
            next_class.record_members == 0 ||
            (!bucket_classes.empty() &&
             !mask_less(bucket_classes.back().representative,
                        next_class.representative))) {
          throw std::runtime_error("An exact class record is inconsistent");
        }
        bucket_classes.push_back(next_class);
        have_class = classes.next(next_class);
      }
      if (bucket_classes.empty() ||
          (have_class && next_class.bucket < repeated_buckets)) {
        throw std::runtime_error("A repeated bucket lacks ordered class records");
      }

      std::vector<std::uint64_t> observed_records(bucket_classes.size());
      std::vector<std::uint64_t> observed_members(bucket_classes.size());
      std::vector<std::uint32_t> representative_class(
          bucket.size(), std::numeric_limits<std::uint32_t>::max());
      for (std::uint32_t class_index = 0;
           class_index < bucket_classes.size(); ++class_index) {
        const auto representative_index =
            bucket_classes[class_index].representative_index;
        if (representative_class[representative_index] !=
            std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error("Two classes use one representative index");
        }
        representative_class[representative_index] = class_index;
      }
      for (std::uint32_t candidate_index = 0;
           candidate_index < bucket.size(); ++candidate_index) {
        AssignmentRecord assignment{};
        if (!assignments.next(assignment) ||
            assignment.local_class >= bucket_classes.size()) {
          throw std::runtime_error("An assignment record is missing or invalid");
        }
        const auto expected_class = representative_class[candidate_index];
        if (expected_class != std::numeric_limits<std::uint32_t>::max() &&
            (assignment.local_class != expected_class || !identity(assignment))) {
          throw std::runtime_error("A class representative lacks its identity witness");
        }
        ++observed_records[assignment.local_class];
        checked_add(observed_members[assignment.local_class],
                    bucket[candidate_index].members);
        ++assignment_records;
      }
      for (std::uint32_t class_index = 0;
           class_index < bucket_classes.size(); ++class_index) {
        if (observed_records[class_index] !=
                bucket_classes[class_index].record_members ||
            observed_members[class_index] == 0) {
          throw std::runtime_error("Weighted class accounting did not close");
        }
        write_record(output,
                     WeightedRecord{signature,
                                    bucket_classes[class_index].representative,
                                    observed_members[class_index]});
        checked_add(output_members, observed_members[class_index]);
        ++exact_classes;
      }
      split_buckets += bucket_classes.size() > 1;
      repeated_records += bucket.size();
      ++repeated_buckets;
      have_repeated = read_fixed(repeated, next_repeated,
                                 "The repeated-bucket index");
    }

    AssignmentRecord extra_assignment{};
    output.close();
    if (!output || have_repeated || have_class ||
        assignments.next(extra_assignment) ||
        source_members != output_members ||
        repeated_records + singleton_classes != source_records ||
        assignment_records != repeated_records ||
        fs::file_size(temporary) !=
            (singleton_classes + exact_classes) * sizeof(WeightedRecord)) {
      throw std::runtime_error("Global weighted representative accounting failed");
    }
    fs::rename(temporary, arguments.output);
    temporary.clear();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"complete_length54_m09_global_weighted_representatives_v1\",\n"
        << "  \"source_record_count\": " << source_records << ",\n"
        << "  \"source_member_count\": " << source_members << ",\n"
        << "  \"signature_bucket_count\": " << signature_buckets << ",\n"
        << "  \"repeated_bucket_count\": " << repeated_buckets << ",\n"
        << "  \"repeated_record_count\": " << repeated_records << ",\n"
        << "  \"singleton_class_count\": " << singleton_classes << ",\n"
        << "  \"repeated_exact_class_count\": " << exact_classes << ",\n"
        << "  \"representative_count\": "
        << (singleton_classes + exact_classes) << ",\n"
        << "  \"output_member_count\": " << output_members << ",\n"
        << "  \"split_bucket_count\": " << split_buckets << ",\n"
        << "  \"assignment_record_count\": " << assignment_records << ",\n"
        << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(arguments.output)
        << ",\n"
        << "  \"assignment_input_count\": " << assignment_paths.size() << ",\n"
        << "  \"class_input_count\": " << class_paths.size() << ",\n"
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
