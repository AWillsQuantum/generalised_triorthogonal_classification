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

constexpr int kDimension = 14;
constexpr int kSupportSize = 54;

struct Mask {
  std::array<std::uint64_t, 256> words{};
};

struct LedgerRecord {
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

struct AssignmentRecord {
  std::uint32_t local_class_index = 0;
  std::uint16_t translation = 0;
  std::array<std::uint16_t, kDimension> standard_basis_images{};
};

struct ClassRecord {
  std::uint64_t bucket_index = 0;
  std::uint32_t local_class_index = 0;
  std::uint32_t representative_candidate_index = 0;
  Mask representative_mask{};
  std::uint64_t member_count = 0;
};

struct NegativeRecord {
  std::uint64_t bucket_index = 0;
  std::uint32_t candidate_index = 0;
  std::uint32_t local_class_index = 0;
};

static_assert(sizeof(LedgerRecord) == 2080);
static_assert(sizeof(BucketRecord) == 4144);
static_assert(sizeof(AssignmentRecord) == 36);
static_assert(sizeof(ClassRecord) == 2072);
static_assert(sizeof(NegativeRecord) == 16);

struct Arguments {
  fs::path ledger;
  fs::path bucket_index;
  std::uint64_t first_bucket = 0;
  std::uint64_t last_bucket = 0;
  bool saw_first = false;
  bool saw_last = false;
  fs::path assignments;
  fs::path classes;
  fs::path negatives;
};

bool mask_equal(const Mask& left, const Mask& right) {
  return left.words == right.words;
}

bool point_is_set(const Mask& mask, int point) {
  return ((mask.words[point / 64] >> (point % 64)) & 1U) != 0;
}

void set_point(Mask& mask, int point) {
  mask.words[point / 64] |= std::uint64_t{1} << (point % 64);
}

std::array<std::uint16_t, kSupportSize> extract_points(const Mask& mask) {
  std::array<std::uint16_t, kSupportSize> result{};
  int count = 0;
  for (int point = 0; point < (1 << kDimension); ++point) {
    if (!point_is_set(mask, point)) continue;
    if (count == kSupportSize) {
      throw std::runtime_error("A quotient mask has weight greater than 54");
    }
    result[count++] = point;
  }
  if (count != kSupportSize) {
    throw std::runtime_error("A quotient mask does not have weight 54");
  }
  return result;
}

int gf2_rank(const std::array<std::uint16_t, kDimension>& vectors) {
  std::array<std::uint16_t, kDimension> basis{};
  int rank = 0;
  for (auto value : vectors) {
    for (int pivot = kDimension - 1; pivot >= 0; --pivot) {
      if (((value >> pivot) & 1U) == 0) continue;
      if (basis[pivot]) {
        value ^= basis[pivot];
      } else {
        basis[pivot] = value;
        ++rank;
        break;
      }
    }
  }
  return rank;
}

std::uint16_t apply_linear(
    std::uint16_t point,
    const std::array<std::uint16_t, kDimension>& images) {
  std::uint16_t result = 0;
  while (point) {
    const int bit = __builtin_ctz(point);
    result ^= images[bit];
    point &= point - 1;
  }
  return result;
}

bool witness_maps_support(
    const std::array<std::uint16_t, kSupportSize>& representative,
    const AssignmentRecord& assignment, const Mask& target) {
  if (gf2_rank(assignment.standard_basis_images) != kDimension) return false;
  Mask image;
  for (const auto point : representative) {
    set_point(image, assignment.translation ^
                         apply_linear(point, assignment.standard_basis_images));
  }
  return mask_equal(image, target);
}

template <typename T>
bool read_record(std::ifstream& input, T& record) {
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == 0 && input.eof()) return false;
  if (input.gcount() != static_cast<std::streamsize>(sizeof(record))) {
    throw std::runtime_error("A quotient verifier input is truncated");
  }
  return true;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string option = argv[index];
    const std::string value = argv[++index];
    if (option == "--ledger") result.ledger = value;
    else if (option == "--bucket-index") result.bucket_index = value;
    else if (option == "--first-bucket") {
      result.first_bucket = std::stoull(value);
      result.saw_first = true;
    } else if (option == "--last-bucket") {
      result.last_bucket = std::stoull(value);
      result.saw_last = true;
    } else if (option == "--assignments") result.assignments = value;
    else if (option == "--classes") result.classes = value;
    else if (option == "--negatives") result.negatives = value;
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.ledger.empty() || result.bucket_index.empty() ||
      result.assignments.empty() || result.classes.empty() ||
      result.negatives.empty() || !result.saw_first || !result.saw_last ||
      result.first_bucket >= result.last_bucket) {
    throw std::runtime_error("Missing or invalid verifier arguments");
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
        !fs::is_regular_file(arguments.assignments) ||
        !fs::is_regular_file(arguments.classes) ||
        !fs::is_regular_file(arguments.negatives) ||
        fs::file_size(arguments.ledger) % sizeof(LedgerRecord) ||
        fs::file_size(arguments.bucket_index) % sizeof(BucketRecord) ||
        fs::file_size(arguments.assignments) % sizeof(AssignmentRecord) ||
        fs::file_size(arguments.classes) % sizeof(ClassRecord) ||
        fs::file_size(arguments.negatives) % sizeof(NegativeRecord)) {
      throw std::runtime_error("A quotient verifier input has a bad size");
    }

    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream buckets(arguments.bucket_index, std::ios::binary);
    std::ifstream assignments(arguments.assignments, std::ios::binary);
    std::ifstream classes(arguments.classes, std::ios::binary);
    std::ifstream negatives(arguments.negatives, std::ios::binary);
    if (!ledger || !buckets || !assignments || !classes || !negatives) {
      throw std::runtime_error("Could not open a quotient verifier stream");
    }
    buckets.seekg(static_cast<std::streamoff>(
        arguments.first_bucket * sizeof(BucketRecord)));

    std::vector<ClassRecord> class_records;
    ClassRecord class_record;
    while (read_record(classes, class_record)) class_records.push_back(class_record);
    std::size_t class_cursor = 0;
    std::uint64_t candidate_count = 0;
    std::uint64_t verified_witness_count = 0;
    std::uint64_t negative_count = 0;
    std::uint64_t split_bucket_count = 0;

    for (std::uint64_t bucket_index = arguments.first_bucket;
         bucket_index < arguments.last_bucket; ++bucket_index) {
      BucketRecord bucket;
      if (!read_record(buckets, bucket) || bucket.candidate_count == 0) {
        throw std::runtime_error("The quotient bucket range is truncated");
      }
      std::vector<ClassRecord> local_classes;
      while (class_cursor < class_records.size() &&
             class_records[class_cursor].bucket_index == bucket_index) {
        const auto& current = class_records[class_cursor++];
        if (current.local_class_index != local_classes.size() ||
            current.representative_candidate_index >= bucket.candidate_count ||
            current.member_count == 0) {
          throw std::runtime_error("A quotient class record is invalid");
        }
        local_classes.push_back(current);
      }
      if (local_classes.empty()) {
        throw std::runtime_error("A quotient bucket has no class record");
      }
      if (local_classes.size() > 1) ++split_bucket_count;
      std::vector<std::array<std::uint16_t, kSupportSize>> representatives;
      std::vector<std::uint64_t> member_counts(local_classes.size(), 0);
      std::vector<bool> representative_seen(local_classes.size(), false);
      representatives.reserve(local_classes.size());
      for (const auto& current : local_classes) {
        representatives.push_back(extract_points(current.representative_mask));
      }

      ledger.seekg(static_cast<std::streamoff>(
          bucket.ledger_offset_records * sizeof(LedgerRecord)));
      for (std::uint32_t candidate_index = 0;
           candidate_index < bucket.candidate_count; ++candidate_index) {
        LedgerRecord candidate;
        AssignmentRecord assignment;
        if (!read_record(ledger, candidate) || !read_record(assignments, assignment) ||
            candidate.signature != bucket.signature ||
            assignment.local_class_index >= local_classes.size() ||
            !witness_maps_support(
                representatives[assignment.local_class_index], assignment,
                candidate.mask)) {
          throw std::runtime_error("An affine assignment witness is invalid");
        }
        const auto class_index = assignment.local_class_index;
        if (local_classes[class_index].representative_candidate_index ==
            candidate_index) {
          if (!mask_equal(local_classes[class_index].representative_mask,
                          candidate.mask)) {
            throw std::runtime_error("A class representative mask is inconsistent");
          }
          representative_seen[class_index] = true;
        }
        ++member_counts[class_index];
        ++verified_witness_count;
        for (std::uint32_t rejected_class = 0; rejected_class < class_index;
             ++rejected_class) {
          NegativeRecord negative;
          if (!read_record(negatives, negative) ||
              negative.bucket_index != bucket_index ||
              negative.candidate_index != candidate_index ||
              negative.local_class_index != rejected_class) {
            throw std::runtime_error("The negative-comparison ledger is not exact");
          }
          ++negative_count;
        }
      }
      for (std::size_t index = 0; index < local_classes.size(); ++index) {
        if (!representative_seen[index] ||
            member_counts[index] != local_classes[index].member_count) {
          throw std::runtime_error("A class member count does not close");
        }
      }
      candidate_count += bucket.candidate_count;
    }

    AssignmentRecord extra_assignment;
    NegativeRecord extra_negative;
    if (class_cursor != class_records.size() || read_record(assignments, extra_assignment) ||
        read_record(negatives, extra_negative) ||
        candidate_count * sizeof(AssignmentRecord) !=
            fs::file_size(arguments.assignments) ||
        negative_count * sizeof(NegativeRecord) !=
            fs::file_size(arguments.negatives)) {
      throw std::runtime_error("The quotient verifier accounting did not close");
    }
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    std::cout << "{\n"
              << "  \"status\": \"verified_length54_m14_exact_quotient_batch_v1\",\n"
              << "  \"first_bucket\": " << arguments.first_bucket << ",\n"
              << "  \"last_bucket\": " << arguments.last_bucket << ",\n"
              << "  \"candidate_count\": " << candidate_count << ",\n"
              << "  \"class_count\": " << class_records.size() << ",\n"
              << "  \"negative_record_count\": " << negative_count << ",\n"
              << "  \"split_bucket_count\": " << split_bucket_count << ",\n"
              << "  \"verified_witness_count\": " << verified_witness_count << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "verify_m14_quotient_batch: " << error.what() << '\n';
    return 1;
  }
}
