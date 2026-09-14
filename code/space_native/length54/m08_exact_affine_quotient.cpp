#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr int kDimension = 8;
constexpr int kAmbientSize = 1 << kDimension;
constexpr int kMaximumSupportSize = 54;
constexpr std::size_t kAssignmentBytes = 16;
constexpr std::size_t kOutputBufferBytes = 16 * 1024 * 1024;

using Mask = std::array<std::uint64_t, 4>;
using Signature = std::array<std::uint64_t, 4>;
using Profile = std::array<std::uint8_t, kMaximumSupportSize>;
using Witness = std::array<std::uint8_t, kDimension + 1>;

struct LedgerRecord {
  Signature signature{};
  Mask mask{};
};

static_assert(sizeof(LedgerRecord) == 64);

struct SearchMetrics {
  std::uint64_t target_origin_count = 0;
  std::uint64_t search_node_count = 0;
  std::uint64_t rejected_candidate_image_count = 0;
};

struct PreparedSupport {
  int support_size = 0;
  std::array<std::uint8_t, kMaximumSupportSize> points{};
  std::array<std::uint8_t, kAmbientSize> point_index{};
  std::array<std::uint8_t, kAmbientSize> difference_counts{};
  std::array<Profile, kMaximumSupportSize> profiles{};
  std::array<std::uint8_t, kAmbientSize - 1> difference_multiset{};
  std::array<Profile, kMaximumSupportSize> profile_multiset{};
  std::uint8_t source_origin = 0;
  Mask translated_mask{};
  std::array<std::uint8_t, kDimension> source_basis{};
  std::array<std::array<std::uint8_t, kAmbientSize>, kDimension + 1>
      source_spans{};
};

struct Arguments {
  fs::path input;
  fs::path output;
  std::uint64_t first_record = 0;
  std::uint64_t record_count = 0;
  std::uint64_t expected_buckets = 0;
  int support_size = 0;
  Signature expected_first_signature{};
  Signature expected_last_signature{};
  bool have_first_signature = false;
  bool have_last_signature = false;
};

bool mask_less(const Mask& left, const Mask& right) {
  for (int index = 3; index >= 0; --index) {
    if (left[index] != right[index]) return left[index] < right[index];
  }
  return false;
}

bool record_less(const LedgerRecord& left, const LedgerRecord& right) {
  if (left.signature != right.signature) return left.signature < right.signature;
  return mask_less(left.mask, right.mask);
}

bool mask_contains(const Mask& mask, std::uint8_t point) {
  return ((mask[point >> 6] >> (point & 63)) & 1ULL) != 0;
}

void mask_insert(Mask& mask, std::uint8_t point) {
  mask[point >> 6] |= 1ULL << (point & 63);
}

Signature parse_signature(const std::string& value) {
  if (value.size() != 64) {
    throw std::runtime_error("A signature must contain exactly 64 hex digits");
  }
  Signature result{};
  for (int index = 0; index < 4; ++index) {
    result[index] = std::stoull(value.substr(16 * index, 16), nullptr, 16);
  }
  return result;
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input") {
      result.input = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--first-record") {
      result.first_record = std::stoull(value);
    } else if (option == "--record-count") {
      result.record_count = std::stoull(value);
    } else if (option == "--expected-buckets") {
      result.expected_buckets = std::stoull(value);
    } else if (option == "--support-size") {
      result.support_size = std::stoi(value);
    } else if (option == "--expected-first-signature") {
      result.expected_first_signature = parse_signature(value);
      result.have_first_signature = true;
    } else if (option == "--expected-last-signature") {
      result.expected_last_signature = parse_signature(value);
      result.have_last_signature = true;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.input.empty() || result.output.empty() || result.record_count == 0 ||
      result.expected_buckets == 0 || result.support_size < 1 ||
      result.support_size > kMaximumSupportSize || !result.have_first_signature ||
      !result.have_last_signature) {
    throw std::runtime_error(
        "Required: --input PATH --output PATH --first-record N "
        "--record-count N --expected-buckets N --support-size N "
        "--expected-first-signature HEX --expected-last-signature HEX");
  }
  return result;
}

std::array<std::uint8_t, kMaximumSupportSize> points_from_mask(
    const Mask& mask, int expected_size) {
  std::array<std::uint8_t, kMaximumSupportSize> result{};
  int count = 0;
  for (int word_index = 0; word_index < 4; ++word_index) {
    std::uint64_t word = mask[word_index];
    while (word != 0) {
      const auto bit = static_cast<int>(__builtin_ctzll(word));
      if (count >= expected_size) {
        throw std::runtime_error("A ledger support has excessive weight");
      }
      result[count++] = static_cast<std::uint8_t>(64 * word_index + bit);
      word &= word - 1;
    }
  }
  if (count != expected_size) {
    throw std::runtime_error("A ledger support has the wrong weight");
  }
  return result;
}

bool profile_less(const Profile& left, const Profile& right) {
  return left < right;
}

bool feature_equal(
    const PreparedSupport& support,
    std::uint8_t left_vector,
    std::uint8_t right_vector) {
  const auto left_point = static_cast<std::uint8_t>(
      support.source_origin ^ left_vector);
  const auto right_point = static_cast<std::uint8_t>(
      support.source_origin ^ right_vector);
  return support.difference_counts[left_vector] ==
             support.difference_counts[right_vector] &&
         support.profiles[support.point_index[left_point]] ==
             support.profiles[support.point_index[right_point]];
}

bool insert_independent(
    std::array<std::uint8_t, kDimension>& linear_basis,
    std::uint8_t vector) {
  auto residual = vector;
  for (int bit = kDimension - 1; bit >= 0; --bit) {
    if (((residual >> bit) & 1U) == 0) continue;
    if (linear_basis[bit] != 0) {
      residual ^= linear_basis[bit];
    } else {
      linear_basis[bit] = residual;
      return true;
    }
  }
  return false;
}

PreparedSupport prepare_support(const Mask& mask, int support_size) {
  PreparedSupport result;
  result.support_size = support_size;
  result.points = points_from_mask(mask, support_size);
  result.point_index.fill(0xFF);
  result.difference_counts.fill(0);
  for (int index = 0; index < support_size; ++index) {
    result.point_index[result.points[index]] = static_cast<std::uint8_t>(index);
    result.profiles[index].fill(0);
  }
  for (int left = 0; left < support_size; ++left) {
    for (int right = left + 1; right < support_size; ++right) {
      result.difference_counts[result.points[left] ^ result.points[right]] += 2;
    }
  }
  for (int left = 0; left < support_size; ++left) {
    int cursor = 0;
    for (int right = 0; right < support_size; ++right) {
      if (left == right) continue;
      result.profiles[left][cursor++] =
          result.difference_counts[result.points[left] ^ result.points[right]];
    }
    std::sort(
        result.profiles[left].begin(),
        result.profiles[left].begin() + support_size - 1);
  }
  for (int value = 1; value < kAmbientSize; ++value) {
    result.difference_multiset[value - 1] = result.difference_counts[value];
  }
  std::sort(
      result.difference_multiset.begin(), result.difference_multiset.end());
  for (int index = 0; index < support_size; ++index) {
    result.profile_multiset[index] = result.profiles[index];
  }
  std::sort(
      result.profile_multiset.begin(),
      result.profile_multiset.begin() + support_size,
      profile_less);

  result.source_origin = result.points[0];
  for (int index = 0; index < support_size; ++index) {
    mask_insert(
        result.translated_mask,
        static_cast<std::uint8_t>(result.points[index] ^ result.source_origin));
  }

  std::array<std::uint8_t, kMaximumSupportSize - 1> vectors{};
  std::array<int, kMaximumSupportSize - 1> frequencies{};
  int vector_count = 0;
  for (int index = 0; index < support_size; ++index) {
    const auto vector = static_cast<std::uint8_t>(
        result.points[index] ^ result.source_origin);
    if (vector != 0) vectors[vector_count++] = vector;
  }
  for (int left = 0; left < vector_count; ++left) {
    for (int right = 0; right < vector_count; ++right) {
      frequencies[left] += feature_equal(result, vectors[left], vectors[right]);
    }
  }
  std::array<int, kMaximumSupportSize - 1> order{};
  for (int index = 0; index < vector_count; ++index) order[index] = index;
  std::sort(order.begin(), order.begin() + vector_count, [&](int left, int right) {
    if (frequencies[left] != frequencies[right]) {
      return frequencies[left] < frequencies[right];
    }
    const auto left_vector = vectors[left];
    const auto right_vector = vectors[right];
    const auto left_count = result.difference_counts[left_vector];
    const auto right_count = result.difference_counts[right_vector];
    if (left_count != right_count) return left_count < right_count;
    const auto& left_profile = result.profiles[result.point_index[
        static_cast<std::uint8_t>(result.source_origin ^ left_vector)]];
    const auto& right_profile = result.profiles[result.point_index[
        static_cast<std::uint8_t>(result.source_origin ^ right_vector)]];
    if (left_profile != right_profile) return left_profile < right_profile;
    return left_vector < right_vector;
  });
  std::array<std::uint8_t, kDimension> linear_basis{};
  int basis_size = 0;
  for (int cursor = 0; cursor < vector_count && basis_size < kDimension; ++cursor) {
    const auto vector = vectors[order[cursor]];
    if (insert_independent(linear_basis, vector)) {
      result.source_basis[basis_size++] = vector;
    }
  }
  if (basis_size != kDimension) {
    throw std::runtime_error("A ledger support does not have full affine rank");
  }
  result.source_spans[0][0] = 0;
  for (int depth = 0; depth < kDimension; ++depth) {
    const int prior_size = 1 << depth;
    for (int index = 0; index < prior_size; ++index) {
      const auto prior = result.source_spans[depth][index];
      result.source_spans[depth + 1][index] = prior;
      result.source_spans[depth + 1][prior_size + index] =
          static_cast<std::uint8_t>(result.source_basis[depth] ^ prior);
    }
  }
  return result;
}

bool preliminary_invariants_equal(
    const PreparedSupport& source, const PreparedSupport& target) {
  return source.support_size == target.support_size &&
         source.difference_multiset == target.difference_multiset &&
         std::equal(
             source.profile_multiset.begin(),
             source.profile_multiset.begin() + source.support_size,
             target.profile_multiset.begin());
}

struct TransporterSearch {
  const PreparedSupport& source;
  const PreparedSupport& target;
  std::uint8_t target_origin = 0;
  Mask target_mask{};
  std::array<std::uint8_t, kMaximumSupportSize - 1> target_vectors{};
  int target_vector_count = 0;
  std::array<std::array<std::uint8_t, kAmbientSize>, kDimension + 1>
      target_spans{};
  std::array<std::uint8_t, kDimension> images{};
  SearchMetrics& metrics;

  bool search(int depth) {
    ++metrics.search_node_count;
    if (depth == kDimension) return true;
    const int span_size = 1 << depth;
    std::bitset<kAmbientSize> occupied;
    for (int index = 0; index < span_size; ++index) {
      occupied.set(target_spans[depth][index]);
    }
    const auto source_vector = source.source_basis[depth];
    const auto source_point = static_cast<std::uint8_t>(
        source.source_origin ^ source_vector);
    const auto& source_point_profile =
        source.profiles[source.point_index[source_point]];
    for (int vector_index = 0; vector_index < target_vector_count; ++vector_index) {
      const auto target_vector = target_vectors[vector_index];
      const auto target_point = static_cast<std::uint8_t>(
          target_origin ^ target_vector);
      if (occupied.test(target_vector) ||
          target.difference_counts[target_vector] !=
              source.difference_counts[source_vector] ||
          target.profiles[target.point_index[target_point]] !=
              source_point_profile) {
        continue;
      }
      bool valid = true;
      for (int index = 0; index < span_size; ++index) {
        const auto source_new = static_cast<std::uint8_t>(
            source_vector ^ source.source_spans[depth][index]);
        const auto target_new = static_cast<std::uint8_t>(
            target_vector ^ target_spans[depth][index]);
        const bool source_member = mask_contains(source.translated_mask, source_new);
        const bool target_member = mask_contains(target_mask, target_new);
        if (source_member != target_member ||
            source.difference_counts[source_new] !=
                target.difference_counts[target_new] ||
            (source_member &&
             source.profiles[source.point_index[static_cast<std::uint8_t>(
                 source.source_origin ^ source_new)]] !=
                 target.profiles[target.point_index[static_cast<std::uint8_t>(
                     target_origin ^ target_new)]])) {
          valid = false;
          ++metrics.rejected_candidate_image_count;
          break;
        }
        target_spans[depth + 1][span_size + index] = target_new;
      }
      if (!valid) continue;
      for (int index = 0; index < span_size; ++index) {
        target_spans[depth + 1][index] = target_spans[depth][index];
      }
      images[depth] = target_vector;
      if (search(depth + 1)) return true;
    }
    return false;
  }
};

bool witness_valid(
    const PreparedSupport& source,
    const PreparedSupport& target,
    const Witness& witness) {
  std::array<std::uint8_t, kDimension> linear_basis{};
  for (int index = 0; index < kDimension; ++index) {
    if (!insert_independent(linear_basis, witness[index + 1])) return false;
  }
  Mask image{};
  for (int point_index = 0; point_index < source.support_size; ++point_index) {
    const auto point = source.points[point_index];
    auto transformed = witness[0];
    for (int bit = 0; bit < kDimension; ++bit) {
      if (((point >> bit) & 1U) != 0) transformed ^= witness[bit + 1];
    }
    mask_insert(image, transformed);
  }
  Mask target_mask{};
  for (int index = 0; index < target.support_size; ++index) {
    mask_insert(target_mask, target.points[index]);
  }
  return image == target_mask;
}

bool transporter(
    const PreparedSupport& source,
    const PreparedSupport& target,
    Witness& witness,
    SearchMetrics& metrics) {
  if (!preliminary_invariants_equal(source, target)) return false;
  const auto& source_profile =
      source.profiles[source.point_index[source.source_origin]];
  for (int origin_index = 0; origin_index < target.support_size; ++origin_index) {
    const auto target_origin = target.points[origin_index];
    if (target.profiles[target.point_index[target_origin]] != source_profile) continue;
    ++metrics.target_origin_count;
    TransporterSearch search{source, target, target_origin, {}, {}, 0, {}, {}, metrics};
    for (int index = 0; index < target.support_size; ++index) {
      const auto vector = static_cast<std::uint8_t>(
          target.points[index] ^ target_origin);
      mask_insert(search.target_mask, vector);
      if (vector != 0) search.target_vectors[search.target_vector_count++] = vector;
    }
    std::sort(
        search.target_vectors.begin(),
        search.target_vectors.begin() + search.target_vector_count);
    search.target_spans[0][0] = 0;
    if (!search.search(0)) continue;

    std::array<std::uint8_t, kAmbientSize> source_coordinates{};
    for (int coefficient = 0; coefficient < kAmbientSize; ++coefficient) {
      source_coordinates[source.source_spans[kDimension][coefficient]] =
          static_cast<std::uint8_t>(coefficient);
    }
    const auto linear_source_origin = search.target_spans[kDimension][
        source_coordinates[source.source_origin]];
    witness[0] = static_cast<std::uint8_t>(target_origin ^ linear_source_origin);
    for (int bit = 0; bit < kDimension; ++bit) {
      witness[bit + 1] = search.target_spans[kDimension][
          source_coordinates[1 << bit]];
    }
    if (!witness_valid(source, target, witness)) {
      throw std::runtime_error("A discovered affine witness failed validation");
    }
    return true;
  }
  return false;
}

Witness identity_witness() {
  Witness result{};
  for (int bit = 0; bit < kDimension; ++bit) {
    result[bit + 1] = static_cast<std::uint8_t>(1 << bit);
  }
  return result;
}

void append_assignment(
    std::vector<char>& buffer, std::uint16_t class_index, const Witness& witness) {
  std::array<char, kAssignmentBytes> record{};
  record[0] = static_cast<char>(class_index & 0xFF);
  record[1] = static_cast<char>((class_index >> 8) & 0xFF);
  for (int index = 0; index < kDimension + 1; ++index) {
    record[2 + index] = static_cast<char>(witness[index]);
  }
  buffer.insert(buffer.end(), record.begin(), record.end());
}

void flush_buffer(std::ofstream& output, std::vector<char>& buffer) {
  if (buffer.empty()) return;
  output.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  if (!output) throw std::runtime_error("Could not write quotient assignments");
  buffer.clear();
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.input)) {
      throw std::runtime_error("The signature ledger does not exist");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite quotient assignments");
    }
    const auto available_records = fs::file_size(arguments.input) / sizeof(LedgerRecord);
    if (fs::file_size(arguments.input) % sizeof(LedgerRecord) != 0 ||
        arguments.first_record + arguments.record_count > available_records) {
      throw std::runtime_error("The requested ledger slice is invalid");
    }
    std::ifstream input(arguments.input, std::ios::binary);
    input.seekg(static_cast<std::streamoff>(
        arguments.first_record * sizeof(LedgerRecord)));
    const auto temporary = arguments.output.string() + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!input || !output) throw std::runtime_error("Could not open quotient files");

    std::vector<char> output_buffer;
    output_buffer.reserve(kOutputBufferBytes + kAssignmentBytes);
    std::vector<PreparedSupport> representatives;
    std::map<std::uint64_t, std::uint64_t> class_count_distribution;
    Signature current_signature{};
    Signature first_signature{};
    Signature last_signature{};
    bool have_signature = false;
    LedgerRecord prior{};
    bool have_prior = false;
    std::uint64_t bucket_count = 0;
    std::uint64_t affine_class_count = 0;
    std::uint64_t comparison_count = 0;
    std::uint64_t equivalent_comparisons = 0;
    std::uint64_t inequivalent_comparisons = 0;
    SearchMetrics metrics;

    const auto close_bucket = [&]() {
      if (!have_signature) return;
      ++class_count_distribution[representatives.size()];
      affine_class_count += representatives.size();
    };

    for (std::uint64_t record_index = 0;
         record_index < arguments.record_count; ++record_index) {
      LedgerRecord record;
      input.read(reinterpret_cast<char*>(&record), sizeof(record));
      if (!input) throw std::runtime_error("The signature ledger slice is truncated");
      if (have_prior && !record_less(prior, record)) {
        throw std::runtime_error("The signature ledger slice is not strictly sorted");
      }
      prior = record;
      have_prior = true;
      if (!have_signature || record.signature != current_signature) {
        close_bucket();
        current_signature = record.signature;
        representatives.clear();
        ++bucket_count;
        if (!have_signature) first_signature = record.signature;
        last_signature = record.signature;
        have_signature = true;
      }

      auto source = prepare_support(record.mask, arguments.support_size);
      Witness witness{};
      std::size_t class_index = 0;
      bool equivalent = false;
      for (; class_index < representatives.size(); ++class_index) {
        ++comparison_count;
        if (transporter(source, representatives[class_index], witness, metrics)) {
          ++equivalent_comparisons;
          equivalent = true;
          break;
        }
        ++inequivalent_comparisons;
      }
      if (!equivalent) {
        class_index = representatives.size();
        representatives.push_back(source);
        witness = identity_witness();
      }
      if (class_index > 0xFFFF) {
        throw std::runtime_error("A signature bucket exceeds the assignment index format");
      }
      if (!witness_valid(source, representatives[class_index], witness)) {
        throw std::runtime_error("An emitted assignment witness failed validation");
      }
      append_assignment(
          output_buffer, static_cast<std::uint16_t>(class_index), witness);
      if (output_buffer.size() >= kOutputBufferBytes) {
        flush_buffer(output, output_buffer);
      }
    }
    close_bucket();
    flush_buffer(output, output_buffer);
    output.close();
    if (!output || bucket_count != arguments.expected_buckets ||
        first_signature != arguments.expected_first_signature ||
        last_signature != arguments.expected_last_signature ||
        comparison_count != equivalent_comparisons + inequivalent_comparisons ||
        fs::file_size(temporary) != arguments.record_count * kAssignmentBytes) {
      fs::remove(temporary);
      throw std::runtime_error("The exact quotient accounting did not close");
    }
    fs::rename(temporary, arguments.output);
    const auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    std::cout << "{\n"
              << "  \"status\": \"complete_m08_native_exact_quotient_batch_v1\",\n"
              << "  \"support_size\": " << arguments.support_size << ",\n"
              << "  \"ambient_dimension\": 8,\n"
              << "  \"first_record\": " << arguments.first_record << ",\n"
              << "  \"candidate_count\": " << arguments.record_count << ",\n"
              << "  \"bucket_count\": " << bucket_count << ",\n"
              << "  \"affine_class_count\": " << affine_class_count << ",\n"
              << "  \"exact_transporter_comparison_count\": "
              << comparison_count << ",\n"
              << "  \"equivalent_comparison_count\": "
              << equivalent_comparisons << ",\n"
              << "  \"inequivalent_comparison_count\": "
              << inequivalent_comparisons << ",\n"
              << "  \"target_origin_count\": " << metrics.target_origin_count
              << ",\n"
              << "  \"search_node_count\": " << metrics.search_node_count
              << ",\n"
              << "  \"rejected_candidate_image_count\": "
              << metrics.rejected_candidate_image_count << ",\n"
              << "  \"bucket_affine_class_count_distribution\": {";
    bool first = true;
    for (const auto& [class_count, frequency] : class_count_distribution) {
      if (!first) std::cout << ',';
      std::cout << "\n    \"" << class_count << "\": " << frequency;
      first = false;
    }
    if (!first) std::cout << '\n';
    std::cout << "  },\n"
              << "  \"assignment_record_bytes\": 16,\n"
              << "  \"assignment_bytes\": " << fs::file_size(arguments.output)
              << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
