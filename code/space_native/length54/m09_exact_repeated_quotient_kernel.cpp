#include <algorithm>
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

constexpr int kDimension = 9;
constexpr int kAmbientSize = 1 << kDimension;
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

using Profile = std::array<std::uint8_t, kSupportSize - 1>;

struct SupportData {
  Mask mask{};
  std::array<std::uint16_t, kSupportSize> points{};
  std::array<std::uint16_t, kAmbientSize> differences{};
  std::array<Profile, kSupportSize> profiles{};
  std::array<std::int16_t, kAmbientSize> profile_indices{};
  std::array<std::uint16_t, kAmbientSize - 1> difference_multiset{};
  std::array<Profile, kSupportSize> profile_multiset{};
};

struct PreparedSource {
  SupportData data{};
  std::uint16_t origin = 0;
  std::array<std::uint16_t, kDimension> basis{};
  std::array<std::array<std::uint16_t, kAmbientSize>, kDimension + 1> spans{};
};

struct SearchMetrics {
  std::uint64_t target_origin_count = 0;
  std::uint64_t search_node_count = 0;
  std::uint64_t rejected_candidate_image_count = 0;
};

struct Witness {
  std::uint16_t translation = 0;
  std::array<std::uint16_t, kDimension> standard_basis_images{};
};

struct Arguments {
  fs::path ledger;
  fs::path repeated_bucket_index;
  std::uint64_t first_bucket = 0;
  std::uint64_t last_bucket = 0;
  fs::path assignments;
  fs::path classes;
  fs::path negatives;
};

static_assert(sizeof(Mask) == 64);
static_assert(sizeof(LedgerRecord) == 96);
static_assert(sizeof(RepeatedBucketRecord) == 16);
static_assert(sizeof(AssignmentRecord) == 24);
static_assert(sizeof(ClassRecord) == 88);
static_assert(sizeof(NegativeRecord) == 16);

bool point_is_set(const Mask& mask, int point) {
  return ((mask.words[point / 64] >> (point % 64)) & 1U) != 0;
}

int mask_weight(const Mask& mask) {
  int result = 0;
  for (const auto word : mask.words) result += __builtin_popcountll(word);
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

bool mask_equal(const Mask& left, const Mask& right) {
  return left.words == right.words;
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

bool add_to_basis(std::array<std::uint16_t, kDimension>& rows,
                  std::uint16_t vector) {
  auto value = vector;
  for (int pivot = kDimension - 1; pivot >= 0; --pivot) {
    if (((value >> pivot) & 1U) == 0) continue;
    if (rows[pivot]) {
      value ^= rows[pivot];
    } else {
      rows[pivot] = value;
      return true;
    }
  }
  return false;
}

SupportData prepare_support(const Mask& mask) {
  if (mask_weight(mask) != kSupportSize) {
    throw std::runtime_error("A quotient candidate does not have weight 54");
  }
  SupportData result;
  result.mask = mask;
  result.profile_indices.fill(-1);
  int count = 0;
  for (int point = 0; point < kAmbientSize; ++point) {
    if (point_is_set(mask, point)) result.points[count++] = point;
  }
  if (count != kSupportSize) throw std::runtime_error("Point extraction failed");
  for (int left = 0; left < kSupportSize; ++left) {
    for (int right = left + 1; right < kSupportSize; ++right) {
      result.differences[result.points[left] ^ result.points[right]] += 2;
    }
  }
  for (int point = 1; point < kAmbientSize; ++point) {
    result.difference_multiset[point - 1] = result.differences[point];
  }
  std::sort(result.difference_multiset.begin(), result.difference_multiset.end());
  for (int row = 0; row < kSupportSize; ++row) {
    int offset = 0;
    for (int column = 0; column < kSupportSize; ++column) {
      if (row == column) continue;
      result.profiles[row][offset++] = static_cast<std::uint8_t>(
          result.differences[result.points[row] ^ result.points[column]]);
    }
    std::sort(result.profiles[row].begin(), result.profiles[row].end());
    result.profile_indices[result.points[row]] = row;
    result.profile_multiset[row] = result.profiles[row];
  }
  std::sort(result.profile_multiset.begin(), result.profile_multiset.end());
  return result;
}

const Profile& profile(const SupportData& support, int point) {
  const int index = support.profile_indices[point];
  if (index < 0) throw std::runtime_error("Requested a profile outside the support");
  return support.profiles[static_cast<std::size_t>(index)];
}

PreparedSource prepare_source(const SupportData& data) {
  PreparedSource result;
  result.data = data;
  result.origin = data.points[0];
  struct Feature {
    std::uint16_t vector = 0;
    std::uint16_t difference_count = 0;
    Profile profile{};
    int frequency = 0;
  };
  std::vector<Feature> features;
  features.reserve(kSupportSize - 1);
  for (int index = 1; index < kSupportSize; ++index) {
    const auto vector = static_cast<std::uint16_t>(data.points[index] ^ result.origin);
    features.push_back(
        {vector, data.differences[vector], profile(data, data.points[index]), 0});
  }
  for (auto& feature : features) {
    feature.frequency = static_cast<int>(std::count_if(
        features.begin(), features.end(), [&](const Feature& other) {
          return feature.difference_count == other.difference_count &&
                 feature.profile == other.profile;
        }));
  }
  std::sort(features.begin(), features.end(), [](const Feature& left,
                                                  const Feature& right) {
    if (left.frequency != right.frequency) return left.frequency < right.frequency;
    if (left.difference_count != right.difference_count) {
      return left.difference_count < right.difference_count;
    }
    if (left.profile != right.profile) return left.profile < right.profile;
    return left.vector < right.vector;
  });
  std::array<std::uint16_t, kDimension> elimination{};
  int selected = 0;
  for (const auto& feature : features) {
    if (add_to_basis(elimination, feature.vector)) {
      result.basis[selected++] = feature.vector;
      if (selected == kDimension) break;
    }
  }
  if (selected != kDimension) {
    throw std::runtime_error("A quotient candidate lacks full affine rank nine");
  }
  result.spans[0][0] = 0;
  for (int depth = 0; depth < kDimension; ++depth) {
    const int old_size = 1 << depth;
    for (int index = 0; index < old_size; ++index) {
      result.spans[depth + 1][index] = result.spans[depth][index];
      result.spans[depth + 1][old_size + index] =
          result.basis[depth] ^ result.spans[depth][index];
    }
  }
  return result;
}

std::uint16_t apply_linear(
    std::uint16_t point,
    const std::array<std::uint16_t, kDimension>& basis_images) {
  std::uint16_t result = 0;
  for (int bit = 0; bit < kDimension; ++bit) {
    if ((point >> bit) & 1U) result ^= basis_images[bit];
  }
  return result;
}

bool witness_is_valid(const PreparedSource& source, const SupportData& target,
                      const Witness& witness) {
  if (gf2_rank(witness.standard_basis_images) != kDimension) return false;
  for (const auto point : source.data.points) {
    const int image = witness.translation ^
                      apply_linear(point, witness.standard_basis_images);
    if (!point_is_set(target.mask, image)) return false;
  }
  return true;
}

bool search_basis_images(
    const PreparedSource& source, const SupportData& target,
    std::uint16_t target_origin,
    const std::array<std::uint16_t, kSupportSize - 1>& target_vectors,
    int depth, std::array<std::uint16_t, kAmbientSize>& target_span,
    std::array<std::uint16_t, kDimension>& selected_images,
    SearchMetrics& metrics) {
  ++metrics.search_node_count;
  if (depth == kDimension) return true;
  const int span_size = 1 << depth;
  std::array<bool, kAmbientSize> occupied{};
  for (int index = 0; index < span_size; ++index) {
    occupied[target_span[index]] = true;
  }
  const auto source_vector = source.basis[depth];
  const auto& source_profile =
      profile(source.data, source.origin ^ source_vector);
  for (const auto target_vector : target_vectors) {
    if (occupied[target_vector] ||
        target.differences[target_vector] !=
            source.data.differences[source_vector] ||
        profile(target, target_origin ^ target_vector) != source_profile) {
      continue;
    }
    bool valid = true;
    for (int index = 0; index < span_size; ++index) {
      const int source_new = source_vector ^ source.spans[depth][index];
      const int target_new = target_vector ^ target_span[index];
      const bool source_member =
          point_is_set(source.data.mask, source.origin ^ source_new);
      const bool target_member = point_is_set(target.mask, target_origin ^ target_new);
      if (source_member != target_member ||
          source.data.differences[source_new] != target.differences[target_new] ||
          (source_member &&
           profile(source.data, source.origin ^ source_new) !=
               profile(target, target_origin ^ target_new))) {
        valid = false;
        ++metrics.rejected_candidate_image_count;
        break;
      }
      target_span[span_size + index] = target_new;
    }
    if (!valid) continue;
    selected_images[depth] = target_vector;
    if (search_basis_images(source, target, target_origin, target_vectors,
                            depth + 1, target_span, selected_images, metrics)) {
      return true;
    }
  }
  return false;
}

bool transporter(const PreparedSource& source, const SupportData& target,
                 Witness& witness, SearchMetrics& metrics) {
  if (source.data.difference_multiset != target.difference_multiset ||
      source.data.profile_multiset != target.profile_multiset) {
    return false;
  }
  const auto& source_origin_profile = profile(source.data, source.origin);
  for (const auto target_origin : target.points) {
    if (profile(target, target_origin) != source_origin_profile) continue;
    ++metrics.target_origin_count;
    std::array<std::uint16_t, kSupportSize - 1> target_vectors{};
    int vector_count = 0;
    for (const auto point : target.points) {
      if (point != target_origin) target_vectors[vector_count++] = point ^ target_origin;
    }
    std::sort(target_vectors.begin(), target_vectors.end());
    std::array<std::uint16_t, kAmbientSize> target_span{};
    std::array<std::uint16_t, kDimension> selected_images{};
    target_span[0] = 0;
    if (!search_basis_images(source, target, target_origin, target_vectors, 0,
                             target_span, selected_images, metrics)) {
      continue;
    }
    const auto& source_full_span = source.spans[kDimension];
    std::array<std::uint16_t, kAmbientSize> target_full_span{};
    target_full_span[0] = 0;
    for (int depth = 0; depth < kDimension; ++depth) {
      const int old_size = 1 << depth;
      for (int index = 0; index < old_size; ++index) {
        target_full_span[old_size + index] =
            selected_images[depth] ^ target_full_span[index];
      }
    }
    for (int bit = 0; bit < kDimension; ++bit) {
      const auto location = std::find(source_full_span.begin(),
                                      source_full_span.end(), 1 << bit);
      if (location == source_full_span.end()) {
        throw std::runtime_error("Source basis span omitted a standard basis vector");
      }
      witness.standard_basis_images[bit] =
          target_full_span[location - source_full_span.begin()];
    }
    const auto origin_location = std::find(source_full_span.begin(),
                                           source_full_span.end(), source.origin);
    if (origin_location == source_full_span.end()) {
      throw std::runtime_error("Source basis span omitted its origin vector");
    }
    witness.translation = target_origin ^
        target_full_span[origin_location - source_full_span.begin()];
    if (!witness_is_valid(source, target, witness)) {
      throw std::runtime_error("Basis search produced an invalid affine witness");
    }
    return true;
  }
  return false;
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
    }
    else if (option == "--first-bucket") {
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

template <typename T>
void write_record(std::ofstream& output, const T& record) {
  output.write(reinterpret_cast<const char*>(&record), sizeof(record));
  if (!output) throw std::runtime_error("Could not write quotient output");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.ledger) ||
        !fs::is_regular_file(arguments.repeated_bucket_index) ||
        fs::file_size(arguments.ledger) % sizeof(LedgerRecord) != 0 ||
        fs::file_size(arguments.repeated_bucket_index) %
                sizeof(RepeatedBucketRecord) !=
            0) {
      throw std::runtime_error("A quotient input has the wrong byte length");
    }
    const auto ledger_record_count =
        fs::file_size(arguments.ledger) / sizeof(LedgerRecord);
    const auto bucket_count = fs::file_size(arguments.repeated_bucket_index) /
                              sizeof(RepeatedBucketRecord);
    if (arguments.last_bucket > bucket_count) {
      throw std::runtime_error("The quotient bucket range is outside the index");
    }
    if (fs::exists(arguments.assignments) || fs::exists(arguments.classes) ||
        fs::exists(arguments.negatives)) {
      throw std::runtime_error("Refusing to overwrite quotient output");
    }
    const auto assignments_temp = arguments.assignments.string() + ".tmp";
    const auto classes_temp = arguments.classes.string() + ".tmp";
    const auto negatives_temp = arguments.negatives.string() + ".tmp";
    std::ofstream assignments(assignments_temp, std::ios::binary | std::ios::trunc);
    std::ofstream classes(classes_temp, std::ios::binary | std::ios::trunc);
    std::ofstream negatives(negatives_temp, std::ios::binary | std::ios::trunc);
    std::ifstream ledger(arguments.ledger, std::ios::binary);
    std::ifstream buckets(arguments.repeated_bucket_index, std::ios::binary);
    if (!assignments || !classes || !negatives || !ledger || !buckets) {
      throw std::runtime_error("Could not open a quotient stream");
    }
    buckets.seekg(static_cast<std::streamoff>(
        arguments.first_bucket * sizeof(RepeatedBucketRecord)));

    std::uint64_t total_candidates = 0;
    std::uint64_t total_classes = 0;
    std::uint64_t total_comparisons = 0;
    std::uint64_t equivalent_comparisons = 0;
    std::uint64_t inequivalent_comparisons = 0;
    SearchMetrics total_metrics;
    std::uint64_t largest_bucket = 0;
    std::uint64_t split_buckets = 0;
    std::uint64_t prior_bucket_end = 0;
    bool has_prior_bucket = false;
    for (std::uint64_t bucket_index = arguments.first_bucket;
         bucket_index < arguments.last_bucket; ++bucket_index) {
      RepeatedBucketRecord bucket;
      buckets.read(reinterpret_cast<char*>(&bucket), sizeof(bucket));
      if (!buckets || bucket.candidate_count <= 1 ||
          bucket.ledger_offset_records + bucket.candidate_count >
              ledger_record_count ||
          (has_prior_bucket &&
           bucket.ledger_offset_records < prior_bucket_end)) {
        throw std::runtime_error(
            "The repeated-bucket index is truncated, overlapping, or invalid");
      }
      prior_bucket_end =
          bucket.ledger_offset_records + bucket.candidate_count;
      has_prior_bucket = true;
      const bool has_previous = bucket.ledger_offset_records != 0;
      ledger.seekg(static_cast<std::streamoff>(
          (bucket.ledger_offset_records - static_cast<std::uint64_t>(has_previous)) *
          sizeof(LedgerRecord)));
      LedgerRecord previous;
      if (has_previous) {
        ledger.read(reinterpret_cast<char*>(&previous), sizeof(previous));
        if (!ledger) throw std::runtime_error("A ledger boundary is truncated");
      }
      std::vector<PreparedSource> representatives;
      std::vector<std::uint32_t> representative_indices;
      std::vector<std::uint64_t> member_counts;
      Mask prior_mask{};
      bool has_prior = false;
      std::array<std::uint64_t, 4> signature{};
      for (std::uint32_t candidate_index = 0;
           candidate_index < bucket.candidate_count; ++candidate_index) {
        LedgerRecord candidate;
        ledger.read(reinterpret_cast<char*>(&candidate), sizeof(candidate));
        if (!ledger) {
          throw std::runtime_error("A repeated ledger bucket is truncated");
        }
        if (candidate_index == 0) {
          signature = candidate.signature;
          if (has_previous && previous.signature == signature) {
            throw std::runtime_error(
                "A repeated-bucket offset does not start a signature bucket");
          }
        }
        if (candidate.signature != signature ||
            (has_prior && !mask_less(prior_mask, candidate.mask))) {
          throw std::runtime_error("A ledger bucket failed its index binding");
        }
        prior_mask = candidate.mask;
        has_prior = true;
        const auto target = prepare_support(candidate.mask);
        std::uint32_t assigned_class = 0;
        Witness assigned_witness;
        bool assigned = false;
        for (std::uint32_t class_index = 0;
             class_index < representatives.size(); ++class_index) {
          SearchMetrics metrics;
          Witness witness;
          const bool equivalent = transporter(
              representatives[class_index], target, witness, metrics);
          ++total_comparisons;
          total_metrics.target_origin_count += metrics.target_origin_count;
          total_metrics.search_node_count += metrics.search_node_count;
          total_metrics.rejected_candidate_image_count +=
              metrics.rejected_candidate_image_count;
          if (equivalent) {
            ++equivalent_comparisons;
            assigned_class = class_index;
            assigned_witness = witness;
            assigned = true;
            break;
          }
          ++inequivalent_comparisons;
          write_record(negatives,
                       NegativeRecord{bucket_index, candidate_index, class_index});
        }
        if (!assigned) {
          assigned_class = static_cast<std::uint32_t>(representatives.size());
          assigned_witness.translation = 0;
          for (int bit = 0; bit < kDimension; ++bit) {
            assigned_witness.standard_basis_images[bit] = 1 << bit;
          }
          representatives.push_back(prepare_source(target));
          representative_indices.push_back(candidate_index);
          member_counts.push_back(0);
        }
        ++member_counts[assigned_class];
        AssignmentRecord assignment;
        assignment.local_class_index = assigned_class;
        assignment.translation = assigned_witness.translation;
        assignment.standard_basis_images = assigned_witness.standard_basis_images;
        write_record(assignments, assignment);
      }
      if (prior_bucket_end < ledger_record_count) {
        LedgerRecord next;
        ledger.read(reinterpret_cast<char*>(&next), sizeof(next));
        if (!ledger || next.signature == signature) {
          throw std::runtime_error(
              "A repeated-bucket count does not end a signature bucket");
        }
      }
      if (representatives.size() > 1) ++split_buckets;
      largest_bucket = std::max(largest_bucket, bucket.candidate_count);
      total_candidates += bucket.candidate_count;
      total_classes += representatives.size();
      for (std::uint32_t class_index = 0;
           class_index < representatives.size(); ++class_index) {
        ClassRecord record;
        record.bucket_index = bucket_index;
        record.local_class_index = class_index;
        record.representative_candidate_index = representative_indices[class_index];
        record.representative_mask = representatives[class_index].data.mask;
        record.member_count = member_counts[class_index];
        write_record(classes, record);
      }
    }
    assignments.close();
    classes.close();
    negatives.close();
    if (!assignments || !classes || !negatives ||
        fs::file_size(assignments_temp) !=
            total_candidates * sizeof(AssignmentRecord) ||
        fs::file_size(classes_temp) != total_classes * sizeof(ClassRecord) ||
        fs::file_size(negatives_temp) !=
            inequivalent_comparisons * sizeof(NegativeRecord) ||
        total_comparisons != equivalent_comparisons + inequivalent_comparisons) {
      fs::remove(assignments_temp);
      fs::remove(classes_temp);
      fs::remove(negatives_temp);
      throw std::runtime_error("Exact quotient output accounting did not close");
    }
    fs::rename(assignments_temp, arguments.assignments);
    fs::rename(classes_temp, arguments.classes);
    fs::rename(negatives_temp, arguments.negatives);
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m09_exact_repeated_quotient_range_v1\",\n"
              << "  \"first_repeated_bucket\": " << arguments.first_bucket
              << ",\n"
              << "  \"last_repeated_bucket\": " << arguments.last_bucket
              << ",\n"
              << "  \"repeated_bucket_count\": "
              << (arguments.last_bucket - arguments.first_bucket) << ",\n"
              << "  \"repeated_candidate_count\": " << total_candidates
              << ",\n"
              << "  \"affine_class_count\": " << total_classes << ",\n"
              << "  \"split_bucket_count\": " << split_buckets << ",\n"
              << "  \"largest_bucket\": " << largest_bucket << ",\n"
              << "  \"comparison_count\": " << total_comparisons << ",\n"
              << "  \"equivalent_comparison_count\": "
              << equivalent_comparisons << ",\n"
              << "  \"inequivalent_comparison_count\": "
              << inequivalent_comparisons << ",\n"
              << "  \"target_origin_count\": "
              << total_metrics.target_origin_count << ",\n"
              << "  \"search_node_count\": "
              << total_metrics.search_node_count << ",\n"
              << "  \"rejected_candidate_image_count\": "
              << total_metrics.rejected_candidate_image_count << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << ",\n"
              << "  \"candidates_per_second\": "
              << (elapsed == 0 ? 0 : total_candidates / elapsed) << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
