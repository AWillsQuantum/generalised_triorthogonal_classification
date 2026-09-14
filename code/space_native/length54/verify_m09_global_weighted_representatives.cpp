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
  fs::path source;
  fs::path repeated;
  fs::path assignments;
  fs::path classes;
  fs::path representatives;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[index];
    if (option == "--weighted-source") result.source = value;
    else if (option == "--repeated-bucket-index") result.repeated = value;
    else if (option == "--assignment-input-list") result.assignments = value;
    else if (option == "--class-input-list") result.classes = value;
    else if (option == "--weighted-representatives") result.representatives = value;
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.source.empty() || result.repeated.empty() ||
      result.assignments.empty() || result.classes.empty() ||
      result.representatives.empty()) {
    throw std::runtime_error(
        "Required: --weighted-source PATH --repeated-bucket-index PATH "
        "--assignment-input-list PATH --class-input-list PATH "
        "--weighted-representatives PATH");
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

std::vector<fs::path> paths_from(const fs::path& list,
                                 std::size_t record_bytes) {
  std::ifstream input(list);
  if (!input) throw std::runtime_error("Could not open an input list");
  std::vector<fs::path> result;
  std::string line;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty()) result.emplace_back(line);
  }
  if (result.empty()) throw std::runtime_error("An input list is empty");
  for (const auto& path : result) {
    if (!fs::is_regular_file(path) || fs::file_size(path) % record_bytes) {
      throw std::runtime_error("An input-list shard has the wrong size");
    }
  }
  return result;
}

template <typename Record>
class ConcatenatedStream {
 public:
  explicit ConcatenatedStream(std::vector<fs::path> paths)
      : paths_(std::move(paths)) {}

  bool read(Record& record) {
    for (;;) {
      if (!input_.is_open()) {
        if (position_ == paths_.size()) return false;
        input_.open(paths_[position_++], std::ios::binary);
        if (!input_) throw std::runtime_error("Could not read an input shard");
      }
      input_.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (input_.gcount() == static_cast<std::streamsize>(sizeof(record))) {
        return true;
      }
      if (input_.gcount() != 0 || !input_.eof()) {
        throw std::runtime_error("An input shard is truncated");
      }
      input_.close();
      input_.clear();
    }
  }

 private:
  std::vector<fs::path> paths_;
  std::size_t position_ = 0;
  std::ifstream input_;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool source_less(const WeightedRecord& left, const WeightedRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool identity_witness(const AssignmentRecord& assignment) {
  if (assignment.translation) return false;
  for (int bit = 0; bit < 9; ++bit) {
    if (assignment.basis_images[bit] != (1U << bit)) return false;
  }
  return true;
}

void add_without_overflow(std::uint64_t& total, std::uint64_t value) {
  if (std::numeric_limits<std::uint64_t>::max() - total < value) {
    throw std::runtime_error("Member-count overflow during verification");
  }
  total += value;
}

void require_representative(std::ifstream& output, const Signature& signature,
                            const Mask& mask, std::uint64_t members) {
  WeightedRecord actual{};
  if (!read_fixed(output, actual, "The representative ledger") ||
      actual.signature != signature || actual.mask != mask ||
      actual.members != members) {
    throw std::runtime_error("A weighted representative is incorrect");
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    for (const auto& path : {arguments.source, arguments.repeated,
                             arguments.representatives}) {
      if (!fs::is_regular_file(path)) {
        throw std::runtime_error("A verifier input is missing");
      }
    }
    if (fs::file_size(arguments.source) % sizeof(WeightedRecord) ||
        fs::file_size(arguments.repeated) % sizeof(RepeatedRecord) ||
        fs::file_size(arguments.representatives) % sizeof(WeightedRecord)) {
      throw std::runtime_error("A verifier input has the wrong size");
    }
    ConcatenatedStream<AssignmentRecord> assignments(
        paths_from(arguments.assignments, sizeof(AssignmentRecord)));
    ConcatenatedStream<ClassRecord> classes(
        paths_from(arguments.classes, sizeof(ClassRecord)));
    std::ifstream source(arguments.source, std::ios::binary);
    std::ifstream repeated(arguments.repeated, std::ios::binary);
    std::ifstream output(arguments.representatives, std::ios::binary);
    if (!source || !repeated || !output) {
      throw std::runtime_error("Could not open a verifier stream");
    }

    WeightedRecord item{};
    bool have_item = read_fixed(source, item, "The weighted source");
    WeightedRecord previous{};
    bool have_previous = false;
    RepeatedRecord repeated_item{};
    bool have_repeated = read_fixed(repeated, repeated_item,
                                    "The repeated-bucket index");
    ClassRecord class_item{};
    bool have_class = classes.read(class_item);
    std::uint64_t source_records = 0;
    std::uint64_t source_members = 0;
    std::uint64_t repeated_buckets = 0;
    std::uint64_t repeated_records = 0;
    std::uint64_t singleton_classes = 0;
    std::uint64_t exact_classes = 0;
    std::uint64_t verified_members = 0;
    std::uint64_t verified_assignments = 0;
    std::uint64_t split_buckets = 0;

    while (have_item) {
      const auto offset = source_records;
      const auto signature = item.signature;
      std::vector<WeightedRecord> candidates;
      do {
        if (!item.members ||
            (have_previous && !source_less(previous, item))) {
          throw std::runtime_error("The weighted source ordering is invalid");
        }
        previous = item;
        have_previous = true;
        candidates.push_back(item);
        add_without_overflow(source_members, item.members);
        ++source_records;
        have_item = read_fixed(source, item, "The weighted source");
      } while (have_item && item.signature == signature);

      if (have_repeated && repeated_item.offset < offset) {
        throw std::runtime_error("The repeated index moved backwards");
      }
      if (candidates.size() == 1) {
        if (have_repeated && repeated_item.offset == offset) {
          throw std::runtime_error("The repeated index names a singleton");
        }
        require_representative(output, signature, candidates[0].mask,
                               candidates[0].members);
        add_without_overflow(verified_members, candidates[0].members);
        ++singleton_classes;
        continue;
      }
      if (!have_repeated || repeated_item.offset != offset ||
          repeated_item.count != candidates.size()) {
        throw std::runtime_error("The repeated index misses a source bucket");
      }

      std::vector<ClassRecord> exact;
      while (have_class && class_item.bucket == repeated_buckets) {
        if (class_item.local_class != exact.size() ||
            class_item.representative_index >= candidates.size() ||
            class_item.representative !=
                candidates[class_item.representative_index].mask ||
            !class_item.record_members ||
            (!exact.empty() &&
             !mask_less(exact.back().representative,
                        class_item.representative))) {
          throw std::runtime_error("The class stream is inconsistent");
        }
        exact.push_back(class_item);
        have_class = classes.read(class_item);
      }
      if (exact.empty() || (have_class && class_item.bucket < repeated_buckets)) {
        throw std::runtime_error("The class stream skipped a repeated bucket");
      }
      std::vector<std::uint64_t> record_counts(exact.size());
      std::vector<std::uint64_t> member_counts(exact.size());
      std::vector<std::uint32_t> representative_at(
          candidates.size(), std::numeric_limits<std::uint32_t>::max());
      for (std::uint32_t class_index = 0; class_index < exact.size();
           ++class_index) {
        auto& slot = representative_at[exact[class_index].representative_index];
        if (slot != std::numeric_limits<std::uint32_t>::max()) {
          throw std::runtime_error("Representative indices are not unique");
        }
        slot = class_index;
      }
      for (std::uint32_t candidate_index = 0;
           candidate_index < candidates.size(); ++candidate_index) {
        AssignmentRecord assignment{};
        if (!assignments.read(assignment) ||
            assignment.local_class >= exact.size()) {
          throw std::runtime_error("The assignment stream is invalid");
        }
        const auto representative_class = representative_at[candidate_index];
        if (representative_class !=
                std::numeric_limits<std::uint32_t>::max() &&
            (assignment.local_class != representative_class ||
             !identity_witness(assignment))) {
          throw std::runtime_error("A representative witness is not identity");
        }
        ++record_counts[assignment.local_class];
        add_without_overflow(member_counts[assignment.local_class],
                             candidates[candidate_index].members);
        ++verified_assignments;
      }
      for (std::uint32_t class_index = 0; class_index < exact.size();
           ++class_index) {
        if (record_counts[class_index] != exact[class_index].record_members ||
            !member_counts[class_index]) {
          throw std::runtime_error("A weighted class count is wrong");
        }
        require_representative(output, signature,
                               exact[class_index].representative,
                               member_counts[class_index]);
        add_without_overflow(verified_members, member_counts[class_index]);
        ++exact_classes;
      }
      split_buckets += exact.size() > 1;
      repeated_records += candidates.size();
      ++repeated_buckets;
      have_repeated = read_fixed(repeated, repeated_item,
                                 "The repeated-bucket index");
    }

    AssignmentRecord extra_assignment{};
    WeightedRecord extra_output{};
    if (have_repeated || have_class || assignments.read(extra_assignment) ||
        read_fixed(output, extra_output, "The representative ledger") ||
        source_members != verified_members ||
        repeated_records + singleton_classes != source_records ||
        verified_assignments != repeated_records ||
        fs::file_size(arguments.representatives) !=
            (singleton_classes + exact_classes) * sizeof(WeightedRecord)) {
      throw std::runtime_error("The independent weighted audit did not close");
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_global_weighted_representatives_v1\",\n"
        << "  \"source_record_count\": " << source_records << ",\n"
        << "  \"source_member_count\": " << source_members << ",\n"
        << "  \"repeated_bucket_count\": " << repeated_buckets << ",\n"
        << "  \"repeated_record_count\": " << repeated_records << ",\n"
        << "  \"singleton_class_count\": " << singleton_classes << ",\n"
        << "  \"repeated_exact_class_count\": " << exact_classes << ",\n"
        << "  \"representative_count\": "
        << (singleton_classes + exact_classes) << ",\n"
        << "  \"verified_member_count\": " << verified_members << ",\n"
        << "  \"verified_assignment_count\": " << verified_assignments << ",\n"
        << "  \"split_bucket_count\": " << split_buckets << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
