#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 9;
constexpr int kSupportSize = 54;

struct Mask {
  std::array<std::uint64_t, 8> words{};
};

struct LedgerRecord {
  std::array<std::uint64_t, 4> signature{};
  Mask mask{};
};

struct RepeatedBucketRecord {
  std::uint64_t ledger_offset_records = 0;
  std::uint64_t candidate_count = 0;
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

static_assert(sizeof(LedgerRecord) == 96);
static_assert(sizeof(RepeatedBucketRecord) == 16);
static_assert(sizeof(AssignmentRecord) == 24);
static_assert(sizeof(ClassRecord) == 88);
static_assert(sizeof(NegativeRecord) == 16);

struct Arguments {
  fs::path ledger;
  fs::path repeated_bucket_index;
  std::uint64_t first_bucket = 0;
  std::uint64_t last_bucket = 0;
  fs::path assignments;
  fs::path classes;
  fs::path negatives;
};

bool mask_equal(const Mask& left, const Mask& right) {
  return left.words == right.words;
}

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 7; index >= 0; --index) {
    if (left.words[index] != right.words[index]) {
      return left.words[index] < right.words[index];
    }
  }
  return false;
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
      throw std::runtime_error("A representative has weight greater than 54");
    }
    result[count++] = point;
  }
  if (count != kSupportSize) {
    throw std::runtime_error("A representative does not have weight 54");
  }
  return result;
}

int gf2_rank(const std::array<std::uint16_t, kDimension>& vectors) {
  std::array<std::uint16_t, kDimension> basis{};
  int rank = 0;
  for (auto vector : vectors) {
    auto value = vector;
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
  auto remaining = point;
  while (remaining) {
    const auto bit = static_cast<int>(__builtin_ctz(remaining));
    result ^= images[bit];
    remaining &= remaining - 1;
  }
  return result;
}

bool witness_maps_support(
    const std::array<std::uint16_t, kSupportSize>& representative,
    const AssignmentRecord& assignment, const Mask& target) {
  if (gf2_rank(assignment.standard_basis_images) != kDimension) return false;
  Mask image;
  for (const auto point : representative) {
    const int mapped = assignment.translation ^
        apply_linear(point, assignment.standard_basis_images);
    set_point(image, mapped);
  }
  return mask_equal(image, target);
}

bool identity_witness(const AssignmentRecord& assignment) {
  if (assignment.translation != 0) return false;
  for (int bit = 0; bit < kDimension; ++bit) {
    if (assignment.standard_basis_images[bit] != (1U << bit)) return false;
  }
  return true;
}

template <typename Record>
std::optional<Record> read_optional(std::ifstream& input) {
  Record record;
  input.read(reinterpret_cast<char*>(&record), sizeof(record));
  if (input.gcount() == static_cast<std::streamsize>(sizeof(record))) {
    return record;
  }
  if (input.gcount() != 0) {
    throw std::runtime_error("A verifier input is truncated");
  }
  return std::nullopt;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  bool saw_first = false;
  bool saw_last = false;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--ledger") result.ledger = value;
    else if (option == "--repeated-bucket-index") {
      result.repeated_bucket_index = value;
    } else if (option == "--first-bucket") {
      result.first_bucket = std::stoull(value);
      saw_first = true;
    } else if (option == "--last-bucket") {
      result.last_bucket = std::stoull(value);
      saw_last = true;
    } else if (option == "--assignments") result.assignments = value;
    else if (option == "--classes") result.classes = value;
    else if (option == "--negatives") result.negatives = value;
    else throw std::runtime_error("Unknown option: " + option);
  }
  if (result.ledger.empty() || result.repeated_bucket_index.empty() ||
      result.assignments.empty() || result.classes.empty() ||
      result.negatives.empty() || !saw_first || !saw_last ||
      result.first_bucket >= result.last_bucket) {
    throw std::runtime_error(
        "Required: --ledger PATH --repeated-bucket-index PATH "
        "--first-bucket N --last-bucket N --assignments PATH "
        "--classes PATH --negatives PATH");
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    for (const auto& path : {arguments.ledger, arguments.repeated_bucket_index,
                             arguments.assignments, arguments.classes,
                             arguments.negatives}) {
      if (!fs::is_regular_file(path)) {
        throw std::runtime_error("A verifier input is missing");
      }
    }
    if (fs::file_size(arguments.ledger) % sizeof(LedgerRecord) ||
        fs::file_size(arguments.repeated_bucket_index) %
            sizeof(RepeatedBucketRecord) ||
        fs::file_size(arguments.assignments) % sizeof(AssignmentRecord) ||
        fs::file_size(arguments.classes) % sizeof(ClassRecord) ||
        fs::file_size(arguments.negatives) % sizeof(NegativeRecord)) {
      throw std::runtime_error("A verifier input has the wrong byte length");
    }

    const auto ledger_record_count =
        fs::file_size(arguments.ledger) / sizeof(LedgerRecord);
    const auto repeated_bucket_count =
        fs::file_size(arguments.repeated_bucket_index) /
        sizeof(RepeatedBucketRecord);
    const auto assignment_record_count =
        fs::file_size(arguments.assignments) / sizeof(AssignmentRecord);
    const auto class_record_count =
        fs::file_size(arguments.classes) / sizeof(ClassRecord);
    const auto negative_record_count =
        fs::file_size(arguments.negatives) / sizeof(NegativeRecord);
    if (arguments.last_bucket > repeated_bucket_count) {
      throw std::runtime_error("The verifier bucket range is outside the index");
    }

    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream buckets(arguments.repeated_bucket_index, std::ios::binary);
    std::ifstream assignments(arguments.assignments, std::ios::binary);
    std::ifstream classes(arguments.classes, std::ios::binary);
    std::ifstream negatives(arguments.negatives, std::ios::binary);
    if (!ledger || !buckets || !assignments || !classes || !negatives) {
      throw std::runtime_error("Could not open a verifier input");
    }
    buckets.seekg(static_cast<std::streamoff>(
        arguments.first_bucket * sizeof(RepeatedBucketRecord)));

    auto next_class = read_optional<ClassRecord>(classes);
    auto next_negative = read_optional<NegativeRecord>(negatives);
    std::uint64_t verified_witnesses = 0;
    std::uint64_t verified_negative_structure = 0;
    std::uint64_t observed_classes = 0;
    std::uint64_t split_bucket_count = 0;
    std::uint64_t expected_assignments = 0;
    std::uint64_t previous_bucket_end = 0;
    bool has_previous_bucket = false;

    for (std::uint64_t bucket_index = arguments.first_bucket;
         bucket_index < arguments.last_bucket; ++bucket_index) {
      RepeatedBucketRecord bucket;
      buckets.read(reinterpret_cast<char*>(&bucket), sizeof(bucket));
      if (!buckets || bucket.candidate_count <= 1 ||
          bucket.ledger_offset_records + bucket.candidate_count >
              ledger_record_count ||
          (has_previous_bucket &&
           bucket.ledger_offset_records < previous_bucket_end)) {
        throw std::runtime_error("The repeated-bucket authority is invalid");
      }
      previous_bucket_end =
          bucket.ledger_offset_records + bucket.candidate_count;
      has_previous_bucket = true;
      expected_assignments += bucket.candidate_count;

      std::vector<ClassRecord> bucket_classes;
      while (next_class && next_class->bucket_index == bucket_index) {
        if (next_class->local_class_index != bucket_classes.size() ||
            next_class->representative_candidate_index >=
                bucket.candidate_count ||
            next_class->member_count == 0 ||
            (bucket_classes.empty() &&
             next_class->representative_candidate_index != 0) ||
            (!bucket_classes.empty() &&
             next_class->representative_candidate_index <=
                 bucket_classes.back().representative_candidate_index)) {
          throw std::runtime_error("An exact-class record is inconsistent");
        }
        bucket_classes.push_back(*next_class);
        next_class = read_optional<ClassRecord>(classes);
      }
      if (bucket_classes.empty() ||
          (next_class && next_class->bucket_index < bucket_index)) {
        throw std::runtime_error("Exact-class records do not cover a bucket");
      }
      observed_classes += bucket_classes.size();
      split_bucket_count += bucket_classes.size() > 1;

      std::optional<LedgerRecord> previous;
      if (bucket.ledger_offset_records != 0) {
        ledger.seekg(static_cast<std::streamoff>(
            (bucket.ledger_offset_records - 1) * sizeof(LedgerRecord)));
        LedgerRecord record;
        ledger.read(reinterpret_cast<char*>(&record), sizeof(record));
        if (!ledger) throw std::runtime_error("A ledger boundary is truncated");
        previous = record;
      } else {
        ledger.seekg(0);
      }

      std::vector<LedgerRecord> candidates(bucket.candidate_count);
      std::vector<AssignmentRecord> bucket_assignments(bucket.candidate_count);
      for (std::uint64_t candidate_index = 0;
           candidate_index < bucket.candidate_count; ++candidate_index) {
        ledger.read(reinterpret_cast<char*>(&candidates[candidate_index]),
                    sizeof(LedgerRecord));
        assignments.read(
            reinterpret_cast<char*>(&bucket_assignments[candidate_index]),
            sizeof(AssignmentRecord));
        if (!ledger || !assignments ||
            (candidate_index > 0 &&
             (candidates[candidate_index].signature !=
                  candidates.front().signature ||
              !mask_less(candidates[candidate_index - 1].mask,
                         candidates[candidate_index].mask)))) {
          throw std::runtime_error("A repeated candidate stream is invalid");
        }
      }
      if (previous && previous->signature == candidates.front().signature) {
        throw std::runtime_error("A repeated bucket does not start at a signature boundary");
      }
      if (previous_bucket_end < ledger_record_count) {
        LedgerRecord next;
        ledger.read(reinterpret_cast<char*>(&next), sizeof(next));
        if (!ledger || next.signature == candidates.front().signature) {
          throw std::runtime_error("A repeated bucket does not end at a signature boundary");
        }
      }

      std::vector<std::array<std::uint16_t, kSupportSize>> representatives;
      representatives.reserve(bucket_classes.size());
      for (const auto& exact_class : bucket_classes) {
        const auto representative_index =
            exact_class.representative_candidate_index;
        const auto& representative_assignment =
            bucket_assignments[representative_index];
        if (!mask_equal(exact_class.representative_mask,
                        candidates[representative_index].mask) ||
            representative_assignment.local_class_index !=
                exact_class.local_class_index ||
            !identity_witness(representative_assignment)) {
          throw std::runtime_error("An exact-class representative is inconsistent");
        }
        representatives.push_back(extract_points(exact_class.representative_mask));
      }

      std::vector<std::uint64_t> observed_member_counts(bucket_classes.size());
      for (std::uint32_t candidate_index = 0;
           candidate_index < bucket.candidate_count; ++candidate_index) {
        const auto& assignment = bucket_assignments[candidate_index];
        if (assignment.local_class_index >= representatives.size() ||
            candidate_index <
                bucket_classes[assignment.local_class_index]
                    .representative_candidate_index ||
            !witness_maps_support(representatives[assignment.local_class_index],
                                  assignment,
                                  candidates[candidate_index].mask)) {
          throw std::runtime_error("An affine assignment witness is invalid");
        }
        ++observed_member_counts[assignment.local_class_index];
        for (std::uint32_t prior_class = 0;
             prior_class < assignment.local_class_index; ++prior_class) {
          if (!next_negative || next_negative->bucket_index != bucket_index ||
              next_negative->candidate_index != candidate_index ||
              next_negative->local_class_index != prior_class) {
            throw std::runtime_error("Negative-comparison structure is incomplete");
          }
          ++verified_negative_structure;
          next_negative = read_optional<NegativeRecord>(negatives);
        }
        if (next_negative && next_negative->bucket_index == bucket_index &&
            next_negative->candidate_index == candidate_index) {
          throw std::runtime_error("Negative-comparison structure has extra records");
        }
        ++verified_witnesses;
      }
      for (std::size_t class_index = 0; class_index < bucket_classes.size();
           ++class_index) {
        if (observed_member_counts[class_index] !=
            bucket_classes[class_index].member_count) {
          throw std::runtime_error("Exact-class member accounting is inconsistent");
        }
      }
    }

    AssignmentRecord extra_assignment;
    if (assignments.read(reinterpret_cast<char*>(&extra_assignment),
                         sizeof(extra_assignment)) ||
        next_class || next_negative || expected_assignments != assignment_record_count ||
        verified_witnesses != assignment_record_count ||
        observed_classes != class_record_count ||
        verified_negative_structure != negative_record_count) {
      throw std::runtime_error("Independent verifier stream accounting did not close");
    }

    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout
        << "{\n"
        << "  \"status\": \"verified_length54_m09_repeated_quotient_batch_v1\",\n"
        << "  \"first_repeated_bucket\": " << arguments.first_bucket << ",\n"
        << "  \"last_repeated_bucket\": " << arguments.last_bucket << ",\n"
        << "  \"repeated_bucket_count\": "
        << (arguments.last_bucket - arguments.first_bucket) << ",\n"
        << "  \"repeated_candidate_count\": " << assignment_record_count << ",\n"
        << "  \"class_count\": " << class_record_count << ",\n"
        << "  \"split_bucket_count\": " << split_bucket_count << ",\n"
        << "  \"negative_record_count\": " << negative_record_count << ",\n"
        << "  \"verified_negative_structure_count\": "
        << verified_negative_structure << ",\n"
        << "  \"verified_witness_count\": " << verified_witnesses << ",\n"
        << "  \"elapsed_seconds\": " << elapsed << ",\n"
        << "  \"witnesses_per_second\": "
        << (elapsed == 0 ? 0 : verified_witnesses / elapsed) << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
