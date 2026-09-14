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

using Signature = std::array<std::uint64_t, 4>;
using Mask = std::array<std::uint64_t, 8>;

struct SourceRecord {
  Signature signature{};
  Mask mask{};
};

struct RepeatedRecord {
  std::uint64_t offset = 0;
  std::uint64_t count = 0;
};

struct ExactClassRecord {
  std::uint64_t bucket = 0;
  std::uint32_t local_class = 0;
  std::uint32_t representative_index = 0;
  Mask representative{};
  std::uint64_t members = 0;
};

struct WeightedRecord {
  Signature signature{};
  Mask representative{};
  std::uint64_t members = 0;
};

static_assert(sizeof(SourceRecord) == 96);
static_assert(sizeof(RepeatedRecord) == 16);
static_assert(sizeof(ExactClassRecord) == 88);
static_assert(sizeof(WeightedRecord) == 104);

struct Arguments {
  fs::path source;
  fs::path repeated;
  fs::path classes;
  fs::path weighted;
};

Arguments arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[index];
    if (option == "--ledger") {
      result.source = value;
    } else if (option == "--repeated-bucket-index") {
      result.repeated = value;
    } else if (option == "--class-input-list") {
      result.classes = value;
    } else if (option == "--weighted-ledger") {
      result.weighted = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.source.empty() || result.repeated.empty() ||
      result.classes.empty() || result.weighted.empty()) {
    throw std::runtime_error(
        "Required: --ledger PATH --repeated-bucket-index PATH "
        "--class-input-list PATH --weighted-ledger PATH");
  }
  return result;
}

std::vector<fs::path> class_paths(const fs::path& list) {
  std::ifstream input(list);
  if (!input) throw std::runtime_error("Could not read the class list");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.empty()) throw std::runtime_error("The class list is empty");
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

class ExactClassStream {
 public:
  explicit ExactClassStream(std::vector<fs::path> inputs)
      : inputs_(std::move(inputs)) {
    for (const auto& path : inputs_) {
      if (!fs::is_regular_file(path) ||
          fs::file_size(path) % sizeof(ExactClassRecord)) {
        throw std::runtime_error("An exact-class shard has the wrong size");
      }
    }
  }

  bool next(ExactClassRecord& record) {
    while (true) {
      if (!current_.is_open()) {
        if (position_ == inputs_.size()) return false;
        current_.open(inputs_[position_++], std::ios::binary);
        if (!current_) throw std::runtime_error("Could not open a class shard");
      }
      current_.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (current_.gcount() == static_cast<std::streamsize>(sizeof(record))) {
        return true;
      }
      if (current_.gcount() != 0 || !current_.eof()) {
        throw std::runtime_error("An exact-class shard is truncated");
      }
      current_.close();
      current_.clear();
    }
  }

 private:
  std::vector<fs::path> inputs_;
  std::size_t position_ = 0;
  std::ifstream current_;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool source_less(const SourceRecord& left, const SourceRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

void require_weighted(std::ifstream& input, const Signature& signature,
                      const Mask& mask, std::uint64_t members) {
  WeightedRecord actual;
  if (!read_fixed(input, actual, "The weighted ledger") ||
      actual.signature != signature || actual.representative != mask ||
      actual.members != members) {
    throw std::runtime_error("A weighted representative record is incorrect");
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto paths = arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(paths.source) ||
        fs::file_size(paths.source) % sizeof(SourceRecord) ||
        !fs::is_regular_file(paths.repeated) ||
        fs::file_size(paths.repeated) % sizeof(RepeatedRecord) ||
        !fs::is_regular_file(paths.weighted) ||
        fs::file_size(paths.weighted) % sizeof(WeightedRecord)) {
      throw std::runtime_error("A weighted-ledger audit input has the wrong size");
    }
    const auto inputs = class_paths(paths.classes);
    ExactClassStream exact(inputs);
    std::ifstream source(paths.source, std::ios::binary);
    std::ifstream repeated(paths.repeated, std::ios::binary);
    std::ifstream weighted(paths.weighted, std::ios::binary);
    if (!source || !repeated || !weighted) {
      throw std::runtime_error("Could not open a weighted-ledger audit stream");
    }

    SourceRecord item{};
    SourceRecord prior{};
    bool have_item = read_fixed(source, item, "The source ledger");
    bool have_prior = false;
    RepeatedRecord repeated_item{};
    bool have_repeated = read_fixed(repeated, repeated_item,
                                    "The repeated-bucket index");
    ExactClassRecord exact_item{};
    bool have_exact = exact.next(exact_item);
    std::uint64_t source_records = 0;
    std::uint64_t buckets = 0;
    std::uint64_t repeated_buckets = 0;
    std::uint64_t repeated_candidates = 0;
    std::uint64_t singleton_classes = 0;
    std::uint64_t classes = 0;
    std::uint64_t members = 0;
    std::uint64_t split_buckets = 0;

    while (have_item) {
      const auto offset = source_records;
      const auto signature = item.signature;
      std::vector<SourceRecord> bucket;
      do {
        if (have_prior && !source_less(prior, item)) {
          throw std::runtime_error("The source ledger is not strictly ordered");
        }
        prior = item;
        have_prior = true;
        bucket.push_back(item);
        ++source_records;
        have_item = read_fixed(source, item, "The source ledger");
      } while (have_item && item.signature == signature);
      ++buckets;

      if (bucket.size() == 1) {
        if (have_repeated && repeated_item.offset == offset) {
          throw std::runtime_error("The repeated index contains a singleton");
        }
        require_weighted(weighted, signature, bucket[0].mask, 1);
        ++singleton_classes;
        ++classes;
        ++members;
        continue;
      }

      if (!have_repeated || repeated_item.offset != offset ||
          repeated_item.count != bucket.size()) {
        throw std::runtime_error("A repeated bucket is absent from its index");
      }
      if (!have_exact || exact_item.bucket != repeated_buckets) {
        throw std::runtime_error("A repeated bucket has no exact class");
      }
      std::uint32_t local_class = 0;
      std::uint64_t bucket_members = 0;
      std::uint64_t bucket_classes = 0;
      Mask previous_representative{};
      bool have_previous_representative = false;
      while (have_exact && exact_item.bucket == repeated_buckets) {
        if (exact_item.local_class != local_class ||
            exact_item.representative_index >= bucket.size() ||
            exact_item.representative !=
                bucket[exact_item.representative_index].mask ||
            exact_item.members == 0 ||
            (have_previous_representative &&
             !mask_less(previous_representative, exact_item.representative))) {
          throw std::runtime_error("An exact-class row is invalid");
        }
        require_weighted(weighted, signature, exact_item.representative,
                         exact_item.members);
        previous_representative = exact_item.representative;
        have_previous_representative = true;
        bucket_members += exact_item.members;
        ++bucket_classes;
        ++classes;
        ++local_class;
        have_exact = exact.next(exact_item);
      }
      if (bucket_members != bucket.size()) {
        throw std::runtime_error("Exact-class weights do not partition a bucket");
      }
      split_buckets += bucket_classes > 1;
      members += bucket_members;
      repeated_candidates += bucket.size();
      ++repeated_buckets;
      have_repeated = read_fixed(repeated, repeated_item,
                                 "The repeated-bucket index");
    }

    WeightedRecord trailing;
    if (have_repeated || have_exact ||
        read_fixed(weighted, trailing, "The weighted ledger") ||
        fs::file_size(paths.weighted) != classes * sizeof(WeightedRecord) ||
        members != source_records ||
        repeated_candidates + singleton_classes != source_records) {
      throw std::runtime_error("The independent weighted audit did not close");
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_local_weighted_representatives_v1\",\n"
        << "  \"source_record_count\": " << source_records << ",\n"
        << "  \"signature_bucket_count\": " << buckets << ",\n"
        << "  \"repeated_bucket_count\": " << repeated_buckets << ",\n"
        << "  \"repeated_candidate_count\": " << repeated_candidates << ",\n"
        << "  \"singleton_class_count\": " << singleton_classes << ",\n"
        << "  \"affine_class_count\": " << classes << ",\n"
        << "  \"split_bucket_count\": " << split_buckets << ",\n"
        << "  \"source_member_count\": " << members << ",\n"
        << "  \"weighted_ledger_bytes\": " << fs::file_size(paths.weighted)
        << ",\n"
        << "  \"class_input_count\": " << inputs.size() << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
