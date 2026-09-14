#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Mask {
  std::array<std::uint64_t, 16> words{};
};

struct LedgerRecord {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

struct RepeatedBucketRecord {
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t candidate_count = 0;
};

struct ClassRecord {
  std::uint64_t bucket_index = 0;
  std::uint32_t local_class_index = 0;
  std::uint32_t representative_candidate_index = 0;
  Mask representative_mask{};
  std::uint64_t member_count = 0;
};

struct WeightedRecord {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
  std::uint64_t member_count = 0;
};

static_assert(sizeof(LedgerRecord) == 160);
static_assert(sizeof(RepeatedBucketRecord) == 16);
static_assert(sizeof(ClassRecord) == 152);
static_assert(sizeof(WeightedRecord) == 168);

struct Arguments {
  fs::path ledger;
  fs::path repeated_bucket_index;
  fs::path class_input_list;
  fs::path output;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[++index];
    if (option == "--ledger") {
      result.ledger = value;
    } else if (option == "--repeated-bucket-index") {
      result.repeated_bucket_index = value;
    } else if (option == "--class-input-list") {
      result.class_input_list = value;
    } else if (option == "--output") {
      result.output = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.ledger.empty() || result.repeated_bucket_index.empty() ||
      result.class_input_list.empty() || result.output.empty()) {
    throw std::runtime_error(
        "Required: --ledger PATH --repeated-bucket-index PATH "
        "--class-input-list PATH --output PATH");
  }
  return result;
}

std::vector<fs::path> read_paths(const fs::path& list_path) {
  std::ifstream input(list_path);
  if (!input) throw std::runtime_error("Could not open the class input list");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.empty()) throw std::runtime_error("The class input list is empty");
  for (const auto& path : result) {
    if (!fs::is_regular_file(path) || fs::file_size(path) % sizeof(ClassRecord)) {
      throw std::runtime_error("A class input has the wrong byte length");
    }
  }
  return result;
}

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 15; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
}

bool record_less(const LedgerRecord& left, const LedgerRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool read_ledger(std::ifstream& input, LedgerRecord& record) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("The compact ledger is truncated");
  }
  return true;
}

bool read_repeated(std::ifstream& input, RepeatedBucketRecord& record) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("The repeated-bucket index is truncated");
  }
  return true;
}

class ClassStream {
 public:
  explicit ClassStream(std::vector<fs::path> paths)
      : paths_(std::move(paths)) {}

  bool read(ClassRecord& record) {
    while (true) {
      if (!stream_.is_open()) {
        if (next_path_ == paths_.size()) return false;
        stream_.open(paths_[next_path_++], std::ios::binary);
        if (!stream_) throw std::runtime_error("Could not open a class input");
      }
      stream_.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (stream_.gcount() == static_cast<std::streamsize>(sizeof(record))) {
        return true;
      }
      if (stream_.gcount() != 0 || !stream_.eof()) {
        throw std::runtime_error("A class input is truncated");
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

void write_record(std::ofstream& output, const WeightedRecord& record) {
  output.write(reinterpret_cast<const char*>(&record), sizeof(record));
  if (!output) throw std::runtime_error("Could not write a weighted representative");
}

}  // namespace

int main(int argc, char** argv) {
  fs::path temporary;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.ledger) ||
        fs::file_size(arguments.ledger) % sizeof(LedgerRecord) ||
        !fs::is_regular_file(arguments.repeated_bucket_index) ||
        fs::file_size(arguments.repeated_bucket_index) %
            sizeof(RepeatedBucketRecord)) {
      throw std::runtime_error("A compact-merge input has the wrong byte length");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite a weighted ledger");
    }
    const auto class_paths = read_paths(arguments.class_input_list);
    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream repeated(arguments.repeated_bucket_index, std::ios::binary);
    ClassStream classes(class_paths);
    temporary = arguments.output.string() + ".tmp";
    if (fs::exists(temporary)) {
      throw std::runtime_error("A temporary weighted ledger already exists");
    }
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!ledger || !repeated || !output) {
      throw std::runtime_error("Could not open a weighted-ledger stream");
    }

    RepeatedBucketRecord next_repeated{};
    bool have_repeated = read_repeated(repeated, next_repeated);
    ClassRecord next_class{};
    bool have_class = classes.read(next_class);
    LedgerRecord current{};
    bool have_current = read_ledger(ledger, current);
    LedgerRecord global_prior{};
    bool have_global_prior = false;
    std::uint64_t ledger_offset = 0;
    std::uint64_t signature_bucket_count = 0;
    std::uint64_t repeated_bucket_count = 0;
    std::uint64_t repeated_candidate_count = 0;
    std::uint64_t singleton_count = 0;
    std::uint64_t affine_class_count = 0;
    std::uint64_t split_bucket_count = 0;
    std::uint64_t source_member_count = 0;
    std::uint64_t largest_bucket = 0;

    while (have_current) {
      const auto group_offset = ledger_offset;
      const auto signature = current.signature;
      std::vector<LedgerRecord> group;
      do {
        if (have_global_prior && !record_less(global_prior, current)) {
          throw std::runtime_error("The compact ledger is not strictly sorted");
        }
        global_prior = current;
        have_global_prior = true;
        group.push_back(current);
        ++ledger_offset;
        have_current = read_ledger(ledger, current);
      } while (have_current && current.signature == signature);

      ++signature_bucket_count;
      largest_bucket = std::max<std::uint64_t>(largest_bucket, group.size());
      if (group.size() == 1) {
        if (have_repeated && next_repeated.ledger_offset_records == group_offset) {
          throw std::runtime_error("A singleton appears in the repeated index");
        }
        write_record(output, WeightedRecord{signature, group[0].mask, 1});
        ++singleton_count;
        ++affine_class_count;
        ++source_member_count;
        continue;
      }

      if (!have_repeated ||
          next_repeated.ledger_offset_records != group_offset ||
          next_repeated.candidate_count != group.size()) {
        throw std::runtime_error("The repeated index does not match the ledger");
      }
      const auto bucket_index = repeated_bucket_count;
      if (!have_class || next_class.bucket_index != bucket_index) {
        throw std::runtime_error("A repeated bucket has no exact class record");
      }
      std::uint32_t expected_local_class = 0;
      std::uint64_t bucket_members = 0;
      std::uint64_t bucket_classes = 0;
      Mask prior_representative{};
      bool have_prior_representative = false;
      while (have_class && next_class.bucket_index == bucket_index) {
        if (next_class.local_class_index != expected_local_class ||
            next_class.representative_candidate_index >= group.size() ||
            next_class.representative_mask.words !=
                group[next_class.representative_candidate_index].mask.words ||
            next_class.member_count == 0 ||
            (have_prior_representative &&
             !mask_less(prior_representative, next_class.representative_mask))) {
          throw std::runtime_error("An exact class record is inconsistent");
        }
        write_record(output,
                     WeightedRecord{signature, next_class.representative_mask,
                                    next_class.member_count});
        prior_representative = next_class.representative_mask;
        have_prior_representative = true;
        bucket_members += next_class.member_count;
        ++bucket_classes;
        ++affine_class_count;
        ++expected_local_class;
        have_class = classes.read(next_class);
      }
      if (bucket_members != group.size()) {
        throw std::runtime_error("Exact class member counts do not close a bucket");
      }
      split_bucket_count += bucket_classes > 1;
      source_member_count += bucket_members;
      repeated_candidate_count += group.size();
      ++repeated_bucket_count;
      have_repeated = read_repeated(repeated, next_repeated);
    }

    output.close();
    if (!output || have_repeated || have_class ||
        fs::file_size(temporary) != affine_class_count * sizeof(WeightedRecord) ||
        source_member_count != ledger_offset ||
        repeated_candidate_count + singleton_count != ledger_offset) {
      throw std::runtime_error("Weighted representative accounting did not close");
    }
    fs::rename(temporary, arguments.output);
    temporary.clear();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"complete_length54_m10_local_weighted_representatives_v1\",\n"
        << "  \"ledger_record_count\": " << ledger_offset << ",\n"
        << "  \"signature_bucket_count\": " << signature_bucket_count << ",\n"
        << "  \"repeated_bucket_count\": " << repeated_bucket_count << ",\n"
        << "  \"repeated_candidate_count\": " << repeated_candidate_count << ",\n"
        << "  \"singleton_class_count\": " << singleton_count << ",\n"
        << "  \"affine_class_count\": " << affine_class_count << ",\n"
        << "  \"split_bucket_count\": " << split_bucket_count << ",\n"
        << "  \"source_member_count\": " << source_member_count << ",\n"
        << "  \"largest_signature_bucket\": " << largest_bucket << ",\n"
        << "  \"weighted_record_bytes\": 168,\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(arguments.output)
        << ",\n"
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
