#include "utsp/native_core.hpp"
#include "utsp/known_outputs.hpp"
#include "utsp/marked_code_canonical.hpp"
#include "utsp/origin_orbits.hpp"
#include "utsp/puncture_embedding.hpp"
#include "utsp/protocol_catalogue.hpp"
#include "utsp/sparse_completion.hpp"
#include "utsp/support_automorphisms.hpp"
#include "utsp/tensor_canonical.hpp"
#include "utsp/tensor_decomposition.hpp"
#include "utsp/tensor_extensions.hpp"
#include "utsp/tensor_primitives.hpp"
#include "utsp/tensor_unmarking.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cctype>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct SourceRecord {
  std::string id;
  std::string family;
};

struct PointedRecord {
  std::uint32_t source_index = 0;
  std::uint16_t base_length = 0;
  std::uint8_t base_affine_dimension = 0;
  std::uint8_t ambient_dimension = 0;
  std::uint8_t parity_case = 0;
  std::int32_t origin = -1;
  std::vector<std::uint32_t> points;
  std::uint64_t original_record_index =
      std::numeric_limits<std::uint64_t>::max();
};

struct Manifest {
  std::uint32_t format_version = 1;
  std::uint32_t length_limit = 0;
  bool protocol_length_filter = false;
  bool base_support_only = false;
  std::vector<SourceRecord> sources;
  std::vector<PointedRecord> supports;
  std::uint64_t label_cache_fingerprint = 0;
  bool label_cache_fingerprint_valid = false;
};

void fingerprint_word(std::uint64_t& state, std::uint64_t word);
[[nodiscard]] std::uint64_t support_cache_fingerprint(
    const PointedRecord& record);

[[nodiscard]] double seconds(Clock::time_point start, Clock::time_point end) {
  return std::chrono::duration<double>(end - start).count();
}

[[nodiscard]] std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

template <typename Integer>
[[nodiscard]] Integer read_little(std::istream& input) {
  static_assert(std::is_integral_v<Integer>);
  using Unsigned = std::make_unsigned_t<Integer>;
  Unsigned value = 0;
  for (std::size_t byte = 0; byte < sizeof(Integer); ++byte) {
    const int next = input.get();
    if (next == std::char_traits<char>::eof()) {
      throw std::runtime_error("unexpected end of manifest");
    }
    value |= static_cast<Unsigned>(static_cast<unsigned char>(next)) <<
             (8 * byte);
  }
  return static_cast<Integer>(value);
}

[[nodiscard]] std::string read_string(std::istream& input) {
  const auto length = read_little<std::uint16_t>(input);
  std::string value(length, '\0');
  input.read(value.data(), length);
  if (!input) throw std::runtime_error("truncated manifest string");
  return value;
}

[[nodiscard]] Manifest read_manifest(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open manifest: " + path);
  std::array<char, 8> magic{};
  input.read(magic.data(), magic.size());
  constexpr std::array<char, 8> expected{'U', 'T', 'S', 'P', 'T', 'S', '1', '\0'};
  if (!input || magic != expected) {
    throw std::runtime_error("manifest has the wrong magic header");
  }
  const auto version = read_little<std::uint32_t>(input);
  if (version != 1 && version != 2) {
    throw std::runtime_error("unsupported manifest version");
  }
  const auto source_count = read_little<std::uint32_t>(input);
  const auto support_count = read_little<std::uint64_t>(input);
  Manifest manifest;
  manifest.format_version = version;
  manifest.length_limit = read_little<std::uint32_t>(input);
  const auto filter_mode = read_little<std::uint32_t>(input);
  if (filter_mode > 3) throw std::runtime_error("unknown manifest format flags");
  manifest.protocol_length_filter = (filter_mode & 1U) != 0;
  manifest.base_support_only = (filter_mode & 2U) != 0;
  std::uint64_t label_cache_fingerprint = 1469598103934665603ULL;
  fingerprint_word(label_cache_fingerprint, manifest.format_version);
  fingerprint_word(label_cache_fingerprint, manifest.length_limit);
  fingerprint_word(
      label_cache_fingerprint, manifest.protocol_length_filter);
  fingerprint_word(label_cache_fingerprint, manifest.base_support_only);
  fingerprint_word(label_cache_fingerprint, support_count);

  manifest.sources.reserve(source_count);
  for (std::uint32_t index = 0; index < source_count; ++index) {
    manifest.sources.push_back(SourceRecord{read_string(input), read_string(input)});
  }
  if (support_count > std::numeric_limits<std::size_t>::max()) {
    throw std::runtime_error("manifest is too large for this platform");
  }
  manifest.supports.reserve(static_cast<std::size_t>(support_count));
  for (std::uint64_t index = 0; index < support_count; ++index) {
    PointedRecord record;
    record.source_index = read_little<std::uint32_t>(input);
    record.base_length = read_little<std::uint16_t>(input);
    record.base_affine_dimension = read_little<std::uint8_t>(input);
    record.ambient_dimension = read_little<std::uint8_t>(input);
    record.parity_case = read_little<std::uint8_t>(input);
    const auto protocol_length = read_little<std::uint8_t>(input);
    record.origin = read_little<std::int32_t>(input);
    record.original_record_index =
        version >= 2 ? read_little<std::uint64_t>(input) : index;
    if (record.source_index >= source_count) {
      throw std::runtime_error("manifest source index is out of range");
    }
    record.points.reserve(protocol_length);
    for (std::uint32_t point = 0; point < protocol_length; ++point) {
      record.points.push_back(read_little<std::uint32_t>(input));
    }
    fingerprint_word(
        label_cache_fingerprint, support_cache_fingerprint(record));
    manifest.supports.push_back(std::move(record));
  }
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("manifest has trailing bytes");
  }
  manifest.label_cache_fingerprint = label_cache_fingerprint;
  manifest.label_cache_fingerprint_valid = true;
  return manifest;
}

template <typename Integer>
void write_little(std::ostream& output, Integer value) {
  static_assert(std::is_integral_v<Integer>);
  using Unsigned = std::make_unsigned_t<Integer>;
  const auto encoded = static_cast<Unsigned>(value);
  for (std::size_t byte = 0; byte < sizeof(Integer); ++byte) {
    output.put(static_cast<char>((encoded >> (8 * byte)) & 0xffU));
  }
  if (!output) throw std::runtime_error("failed while writing manifest");
}

void write_string(std::ostream& output, std::string_view value) {
  if (value.size() >= (std::size_t{1} << 16U)) {
    throw std::runtime_error("manifest source string is too long");
  }
  write_little(output, static_cast<std::uint16_t>(value.size()));
  output.write(value.data(), static_cast<std::streamsize>(value.size()));
  if (!output) throw std::runtime_error("failed while writing manifest string");
}

void write_manifest(const Manifest& manifest, const std::string& path) {
  namespace fs = std::filesystem;
  const fs::path output_path(path);
  if (!output_path.parent_path().empty()) {
    fs::create_directories(output_path.parent_path());
  }
  const fs::path temporary_path = output_path.string() + ".tmp";
  {
    std::ofstream output(temporary_path, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "could not open temporary manifest: " + temporary_path.string());
    }
    constexpr std::array<char, 8> magic{
        'U', 'T', 'S', 'P', 'T', 'S', '1', '\0'};
    output.write(magic.data(), magic.size());
    if (manifest.format_version != 1 && manifest.format_version != 2) {
      throw std::runtime_error("unsupported output manifest version");
    }
    write_little(output, manifest.format_version);
    write_little(output, static_cast<std::uint32_t>(manifest.sources.size()));
    write_little(output, static_cast<std::uint64_t>(manifest.supports.size()));
    write_little(output, manifest.length_limit);
    write_little(
        output,
        static_cast<std::uint32_t>(manifest.protocol_length_filter) |
            (static_cast<std::uint32_t>(manifest.base_support_only) << 1U));
    for (const auto& source : manifest.sources) {
      write_string(output, source.id);
      write_string(output, source.family);
    }
    for (const auto& record : manifest.supports) {
      if (record.points.size() >= (std::size_t{1} << 8U)) {
        throw std::runtime_error("protocol support is too long for format v1");
      }
      write_little(output, record.source_index);
      write_little(output, record.base_length);
      write_little(output, record.base_affine_dimension);
      write_little(output, record.ambient_dimension);
      write_little(output, record.parity_case);
      write_little(output, static_cast<std::uint8_t>(record.points.size()));
      write_little(output, record.origin);
      if (manifest.format_version >= 2) {
        if (record.original_record_index ==
            std::numeric_limits<std::uint64_t>::max()) {
          throw std::runtime_error(
              "v2 manifest record lacks an original record index");
        }
        write_little(output, record.original_record_index);
      }
      for (const auto point : record.points) write_little(output, point);
    }
    output.flush();
    if (!output) throw std::runtime_error("failed to flush temporary manifest");
  }
  std::error_code error;
  fs::rename(temporary_path, output_path, error);
  if (error) {
    fs::remove(temporary_path);
    throw std::runtime_error("could not promote reduced manifest: " +
                             error.message());
  }
}

constexpr std::array<char, 8> kLabelCacheMagic{
    'U', 'T', 'S', 'P', 'L', 'C', '1', '\0'};

void fingerprint_word(std::uint64_t& state, std::uint64_t word) {
  constexpr std::uint64_t prime = 1099511628211ULL;
  for (unsigned byte = 0; byte < 8; ++byte) {
    state ^= (word >> (8U * byte)) & 0xffU;
    state *= prime;
  }
}

[[nodiscard]] std::uint64_t support_cache_fingerprint(
    const PointedRecord& record) {
  std::uint64_t state = 1469598103934665603ULL;
  fingerprint_word(state, record.source_index);
  fingerprint_word(state, record.base_length);
  fingerprint_word(state, record.base_affine_dimension);
  fingerprint_word(state, record.ambient_dimension);
  fingerprint_word(state, record.parity_case);
  fingerprint_word(state, static_cast<std::uint32_t>(record.origin));
  fingerprint_word(state, record.original_record_index);
  fingerprint_word(state, record.points.size());
  for (const auto point : record.points) fingerprint_word(state, point);
  return state;
}

[[nodiscard]] std::uint64_t manifest_cache_fingerprint(
    const Manifest& manifest) {
  if (manifest.label_cache_fingerprint_valid) {
    return manifest.label_cache_fingerprint;
  }
  std::uint64_t state = 1469598103934665603ULL;
  fingerprint_word(state, manifest.format_version);
  fingerprint_word(state, manifest.length_limit);
  fingerprint_word(state, manifest.protocol_length_filter);
  fingerprint_word(state, manifest.base_support_only);
  fingerprint_word(state, manifest.supports.size());
  for (const auto& record : manifest.supports) {
    fingerprint_word(state, support_cache_fingerprint(record));
  }
  return state;
}

struct QuotientCacheSelection {
  std::uint32_t minimum_distance = 0;
  std::vector<std::uint64_t> offsets;
  std::vector<utsp::Mask> rows;

  [[nodiscard]] utsp::QuotientBasisCacheView view() const {
    return utsp::QuotientBasisCacheView{offsets, rows};
  }
};

[[nodiscard]] QuotientCacheSelection read_quotient_cache(
    const std::string& path,
    const Manifest& manifest,
    std::uint64_t selected_start,
    std::uint64_t selected_count,
    std::uint32_t expected_minimum_distance) {
  if (selected_start > manifest.supports.size() ||
      selected_count > manifest.supports.size() - selected_start) {
    throw std::invalid_argument("label cache selection is outside manifest");
  }
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("could not open label cache: " + path);
  std::array<char, 8> magic{};
  input.read(magic.data(), magic.size());
  if (!input || magic != kLabelCacheMagic) {
    throw std::runtime_error("label cache has the wrong magic header");
  }
  const auto version = read_little<std::uint32_t>(input);
  if (version != 1) throw std::runtime_error("unsupported label cache version");
  QuotientCacheSelection result;
  result.minimum_distance = read_little<std::uint32_t>(input);
  if (result.minimum_distance != expected_minimum_distance) {
    throw std::runtime_error("label cache minimum distance does not match branch");
  }
  const auto support_count = read_little<std::uint64_t>(input);
  if (support_count != manifest.supports.size()) {
    throw std::runtime_error("label cache support count does not match manifest");
  }
  const auto fingerprint = read_little<std::uint64_t>(input);
  if (fingerprint != manifest_cache_fingerprint(manifest)) {
    throw std::runtime_error("label cache fingerprint does not match manifest");
  }
  result.offsets.push_back(0);
  result.offsets.reserve(static_cast<std::size_t>(selected_count) + 1);
  for (std::uint64_t index = 0; index < support_count; ++index) {
    const auto dimension = read_little<std::uint8_t>(input);
    if (dimension > 64) {
      throw std::runtime_error("label cache quotient dimension exceeds 64");
    }
    const bool selected =
        index >= selected_start && index < selected_start + selected_count;
    for (std::uint32_t row = 0; row < dimension; ++row) {
      const auto mask = read_little<utsp::Mask>(input);
      if (selected) result.rows.push_back(mask);
    }
    if (selected) result.offsets.push_back(result.rows.size());
  }
  if (input.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("label cache has trailing bytes");
  }
  if (result.offsets.size() != selected_count + 1) {
    throw std::logic_error("label cache selection count is inconsistent");
  }
  return result;
}

[[nodiscard]] std::string argument_value(
    int argc, char** argv, std::string_view name, bool required = true) {
  for (int index = 2; index + 1 < argc; ++index) {
    if (argv[index] == name) return argv[index + 1];
  }
  if (required) throw std::invalid_argument("missing argument " + std::string(name));
  return {};
}

[[nodiscard]] bool has_flag(int argc, char** argv, std::string_view name) {
  for (int index = 2; index < argc; ++index) {
    if (argv[index] == name) return true;
  }
  return false;
}

[[nodiscard]] std::uint64_t parse_unsigned(std::string_view text) {
  std::size_t consumed = 0;
  const auto value = std::stoull(std::string(text), &consumed, 0);
  if (consumed != text.size()) throw std::invalid_argument("invalid integer argument");
  return value;
}

[[nodiscard]] std::vector<std::uint32_t> parse_points(std::string_view text) {
  std::vector<std::uint32_t> points;
  std::size_t start = 0;
  while (start <= text.size()) {
    const auto comma = text.find(',', start);
    const auto token = text.substr(
        start, comma == std::string_view::npos ? text.size() - start : comma - start);
    if (token.empty()) throw std::invalid_argument("empty point in point list");
    const auto value = parse_unsigned(token);
    if (value > std::numeric_limits<std::uint32_t>::max()) {
      throw std::invalid_argument("point exceeds uint32_t");
    }
    points.push_back(static_cast<std::uint32_t>(value));
    if (comma == std::string_view::npos) break;
    start = comma + 1;
  }
  return points;
}

[[nodiscard]] std::vector<std::uint64_t> parse_masks(std::string_view text) {
  std::vector<std::uint64_t> masks;
  std::size_t start = 0;
  while (start <= text.size()) {
    const auto comma = text.find(',', start);
    const auto token = text.substr(
        start, comma == std::string_view::npos ? text.size() - start : comma - start);
    if (token.empty()) throw std::invalid_argument("empty mask in mask list");
    masks.push_back(parse_unsigned(token));
    if (comma == std::string_view::npos) break;
    start = comma + 1;
  }
  return masks;
}

[[nodiscard]] std::vector<std::vector<std::uint32_t>> parse_basis_generators(
    std::string_view text, std::uint32_t dimension) {
  std::vector<std::vector<std::uint32_t>> generators;
  if (text.empty()) return generators;
  std::size_t start = 0;
  while (start <= text.size()) {
    const auto separator = text.find(';', start);
    const auto token = text.substr(
        start,
        separator == std::string_view::npos ? text.size() - start
                                            : separator - start);
    auto basis = parse_points(token);
    if (basis.size() != dimension) {
      throw std::invalid_argument("basis generator has the wrong dimension");
    }
    generators.push_back(std::move(basis));
    if (separator == std::string_view::npos) break;
    start = separator + 1;
  }
  return generators;
}

void print_mask_array(std::span<const utsp::Mask> masks) {
  std::cout << '[';
  for (std::size_t index = 0; index < masks.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '"' << utsp::hex_mask(masks[index]) << '"';
  }
  std::cout << ']';
}

void write_json_string(std::ostream& output, std::string_view value) {
  output.put('"');
  constexpr char hex[] = "0123456789abcdef";
  for (const unsigned char character : value) {
    switch (character) {
      case '"': output << "\\\""; break;
      case '\\': output << "\\\\"; break;
      case '\b': output << "\\b"; break;
      case '\f': output << "\\f"; break;
      case '\n': output << "\\n"; break;
      case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if (character < 0x20U) {
          output << "\\u00" << hex[character >> 4U]
                 << hex[character & 0x0fU];
        } else {
          output.put(static_cast<char>(character));
        }
    }
  }
  output.put('"');
}

[[nodiscard]] std::string fixed_hex(std::uint64_t value) {
  std::ostringstream output;
  output << "0x" << std::hex << std::setfill('0') << std::setw(16) << value;
  return output.str();
}

[[nodiscard]] std::span<const utsp::Mask> cached_basis(
    const QuotientCacheSelection& cache,
    std::size_t index) {
  const auto begin = cache.offsets.at(index);
  const auto end = cache.offsets.at(index + 1);
  return std::span<const utsp::Mask>(
      cache.rows.data() + begin,
      static_cast<std::size_t>(end - begin));
}

[[nodiscard]] std::vector<utsp::Mask> stabilizer_rows(
    const PointedRecord& record) {
  if (record.points.empty() || record.points.size() > 64 ||
      record.ambient_dimension >= 32) {
    throw std::invalid_argument("manifest support dimensions are invalid");
  }
  std::vector<utsp::Mask> rows;
  rows.reserve(record.ambient_dimension);
  for (std::uint32_t coordinate = 0;
       coordinate < record.ambient_dimension;
       ++coordinate) {
    utsp::Mask row = 0;
    for (std::size_t index = 0; index < record.points.size(); ++index) {
      if ((record.points[index] >> coordinate) & 1U) {
        row |= utsp::Mask{1} << index;
      }
    }
    rows.push_back(row);
  }
  return rows;
}

[[nodiscard]] bool quotient_spans_match(
    const PointedRecord& record,
    std::span<const utsp::Mask> reference,
    std::span<const utsp::Mask> candidate) {
  if (reference.size() != candidate.size()) return false;
  const auto width_mask = record.points.size() == 64
                              ? ~utsp::Mask{0}
                              : (utsp::Mask{1} << record.points.size()) - 1;
  for (const auto row : reference) {
    if (row & ~width_mask) return false;
  }
  for (const auto row : candidate) {
    if (row & ~width_mask) return false;
  }
  const auto stabilizer = stabilizer_rows(record);
  std::vector<utsp::Mask> reference_span = stabilizer;
  reference_span.insert(
      reference_span.end(), reference.begin(), reference.end());
  std::vector<utsp::Mask> candidate_span = stabilizer;
  candidate_span.insert(
      candidate_span.end(), candidate.begin(), candidate.end());
  return utsp::rref_basis(reference_span) ==
         utsp::rref_basis(candidate_span);
}

int differential_label_space_caches(int argc, char** argv) {
  const auto started = Clock::now();
  const auto input_path = argument_value(argc, argv, "--input");
  const auto reference_path =
      argument_value(argc, argv, "--reference-cache");
  const auto candidate_path =
      argument_value(argc, argv, "--candidate-cache");
  const auto minimum_distance = static_cast<std::uint32_t>(parse_unsigned(
      argument_value(argc, argv, "--minimum-distance")));
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto requested_workers = workers_text.empty()
                                     ? std::max(1U, std::thread::hardware_concurrency())
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  if (requested_workers == 0) {
    throw std::invalid_argument("cache differential worker count must be positive");
  }

  const auto manifest = read_manifest(input_path);
  const auto reference = read_quotient_cache(
      reference_path,
      manifest,
      0,
      manifest.supports.size(),
      minimum_distance);
  const auto candidate = read_quotient_cache(
      candidate_path,
      manifest,
      0,
      manifest.supports.size(),
      minimum_distance);
  const auto worker_count = std::min<std::size_t>(
      requested_workers, std::max<std::size_t>(1, manifest.supports.size()));
  std::atomic<std::size_t> next{0};
  std::exception_ptr worker_error;
  std::mutex error_mutex;
  std::vector<std::array<std::uint64_t, 65>> histograms(worker_count);
  std::vector<std::uint64_t> identical_counts(worker_count, 0);
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    workers.emplace_back([&, worker] {
      try {
        while (true) {
          const auto index = next.fetch_add(1);
          if (index >= manifest.supports.size()) return;
          const auto reference_basis = cached_basis(reference, index);
          const auto candidate_basis = cached_basis(candidate, index);
          if (reference_basis.size() >= histograms[worker].size()) {
            throw std::runtime_error(
                "cache quotient dimension exceeds histogram range");
          }
          if (!quotient_spans_match(
                  manifest.supports[index],
                  reference_basis,
                  candidate_basis)) {
            throw std::runtime_error(
                "cache quotient spans differ at support " +
                std::to_string(index));
          }
          ++histograms[worker][reference_basis.size()];
          if (std::equal(
                  reference_basis.begin(),
                  reference_basis.end(),
                  candidate_basis.begin(),
                  candidate_basis.end())) {
            ++identical_counts[worker];
          }
        }
      } catch (...) {
        const std::lock_guard lock(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
        next.store(manifest.supports.size());
      }
    });
  }
  for (auto& worker : workers) worker.join();
  if (worker_error) std::rethrow_exception(worker_error);

  std::array<std::uint64_t, 65> histogram{};
  std::uint64_t byte_identical_bases = 0;
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    byte_identical_bases += identical_counts[worker];
    for (std::size_t dimension = 0; dimension < histogram.size(); ++dimension) {
      histogram[dimension] += histograms[worker][dimension];
    }
  }
  std::cout << "{\"schema\":\"utsp-native-label-cache-differential-v1\"," 
            << "\"status\":\"pass\",\"minimum_distance\":"
            << minimum_distance << ",\"supports_compared\":"
            << manifest.supports.size() << ",\"semantic_span_matches\":"
            << manifest.supports.size() << ",\"byte_identical_bases\":"
            << byte_identical_bases << ",\"workers\":" << worker_count
            << ",\"quotient_dimension_histogram\":[";
  bool first = true;
  for (std::size_t dimension = 0; dimension < histogram.size(); ++dimension) {
    if (histogram[dimension] == 0) continue;
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << dimension << ',' << histogram[dimension] << ']';
  }
  std::cout << "],\"elapsed_seconds\":"
            << seconds(started, Clock::now()) << "}\n";
  return 0;
}

[[nodiscard]] std::string row_bits(utsp::Mask row, std::size_t length) {
  std::string result(length, '0');
  for (std::size_t index = 0; index < length; ++index) {
    if (((row >> index) & 1U) != 0) result[index] = '1';
  }
  return result;
}

int self_test() {
  utsp::LabelSpace zero_forms;
  zero_forms.quotient_basis = {1, 2, 4, 8};
  const std::array<std::uint64_t, 4> expected{15, 35, 15, 1};
  bool okay = true;
  for (std::uint32_t q = 1; q <= 4; ++q) {
    okay &= utsp::enumerate_totally_isotropic_subspaces(zero_forms, q) ==
            expected[q - 1];
    okay &= utsp::count_totally_isotropic_subspaces_parallel(
                zero_forms, q, 3) == expected[q - 1];
  }

  utsp::LabelSpace symplectic = zero_forms;
  symplectic.bilinear_form_rows = {{0b0010, 0b0001, 0b1000, 0b0100}};
  okay &= utsp::enumerate_totally_isotropic_subspaces(symplectic, 1) == 15;
  okay &= utsp::enumerate_totally_isotropic_subspaces(symplectic, 2) == 15;
  okay &= utsp::enumerate_totally_isotropic_subspaces(symplectic, 3) == 0;
  okay &= utsp::count_totally_isotropic_subspaces_parallel(
              symplectic, 1, 3) == 15;
  okay &= utsp::count_totally_isotropic_subspaces_parallel(
              symplectic, 2, 3) == 15;
  okay &= utsp::count_totally_isotropic_subspaces_parallel(
              symplectic, 3, 3) == 0;

  utsp::LabelSpace pencil = zero_forms;
  pencil.bilinear_form_rows = {
      {0b0010, 0b0001, 0, 0},
      {0, 0, 0b1000, 0b0100},
  };
  okay &= utsp::alternating_form_span_allows_isotropic_dimension(pencil, 2);
  okay &= !utsp::alternating_form_span_allows_isotropic_dimension(pencil, 3);

  std::vector<std::uint32_t> hard_support(48);
  for (std::uint32_t point = 0; point < hard_support.size(); ++point) {
    hard_support[point] = point;
  }
  const auto hard_label_space =
      utsp::build_logical_label_space(hard_support, 6, false);
  const auto hard_distance_four =
      utsp::build_logical_label_space(hard_support, 6, false, 4);
  const auto hard_distance_four_rebuilt =
      utsp::build_logical_label_space_from_quotient_basis(
          hard_support,
          6,
          hard_distance_four.quotient_basis,
          false);
  const auto hard_distance_four_derived =
      utsp::restrict_logical_label_space_minimum_distance(
          hard_support,
          6,
          hard_label_space.quotient_basis,
          4,
          false);
  okay &= hard_distance_four_rebuilt.points == hard_distance_four.points;
  okay &= hard_distance_four_rebuilt.coordinate_masks ==
          hard_distance_four.coordinate_masks;
  okay &= hard_distance_four_rebuilt.stabilizer_basis ==
          hard_distance_four.stabilizer_basis;
  okay &= hard_distance_four_rebuilt.quotient_basis ==
          hard_distance_four.quotient_basis;
  okay &= hard_distance_four_rebuilt.bilinear_form_rows ==
          hard_distance_four.bilinear_form_rows;
  auto hard_distance_four_span = hard_distance_four.stabilizer_basis;
  hard_distance_four_span.insert(
      hard_distance_four_span.end(),
      hard_distance_four.quotient_basis.begin(),
      hard_distance_four.quotient_basis.end());
  auto hard_distance_four_derived_span =
      hard_distance_four_derived.stabilizer_basis;
  hard_distance_four_derived_span.insert(
      hard_distance_four_derived_span.end(),
      hard_distance_four_derived.quotient_basis.begin(),
      hard_distance_four_derived.quotient_basis.end());
  okay &= utsp::rref_basis(hard_distance_four_span) ==
          utsp::rref_basis(hard_distance_four_derived_span);
  okay &= hard_distance_four.quotient_dimension() ==
          hard_distance_four_derived.quotient_dimension();
  const auto hard_automorphisms =
      utsp::elementary_support_automorphisms(hard_label_space);
  okay &= hard_label_space.quotient_dimension() == 21;
  okay &= hard_automorphisms.elementary_generator_count == 27;
  okay &= hard_automorphisms.generators.size() == 6;
  okay &= hard_automorphisms.subgroup_order == 10'321'920;
  for (const auto& generator : hard_automorphisms.generators) {
    okay &= generator.ambient_basis_images.size() == 6;
    okay &= generator.quotient_basis_images.size() == 21;
    okay &= utsp::gf2_rank(generator.quotient_basis_images) == 21;
  }

  if (utsp::marked_code_canonicalizer_available()) {
    utsp::LabelSpace diagonal_q8;
    diagonal_q8.points.resize(8);
    for (std::uint32_t index = 0; index < 8; ++index) {
      diagonal_q8.points[index] = index;
      diagonal_q8.quotient_basis.push_back(utsp::Mask{1} << index);
    }
    const auto diagonal_result =
        utsp::search_q8_isotropic_subspaces_by_tensor_radical(
            diagonal_q8, 2);
    okay &= diagonal_result.has_isotropic_subspace;
    okay &= diagonal_result.sector_witness_found[8];

    const std::array<utsp::ProtocolSupport, 1> marked_distance_supports{
        utsp::ProtocolSupport{
            846,
            7,
            "self_test_marked_distance",
            "self_test",
            32,
            6,
            0,
            6,
            {0,  1,  8,  11, 13, 15, 18, 19, 20, 21, 22,
             23, 25, 26, 28, 30, 32, 33, 40, 43, 45, 47,
             50, 51, 52, 53, 54, 55, 57, 58, 60, 62},
        }};
    const auto classify_marked_distance = [&](std::uint32_t workers) {
      return utsp::classify_protocol_frontier(
          marked_distance_supports,
          1,
          3,
          24,
          utsp::ProtocolEnumerationMode::marked_code_orbits,
          14,
          workers,
          1,
          0,
          0,
          std::numeric_limits<std::uint32_t>::max(),
          true);
    };
    const auto marked_serial = classify_marked_distance(1);
    const auto marked_parallel = classify_marked_distance(4);
    okay &= marked_serial.marked_code_orbits == 2;
    okay &= marked_serial.nondegenerate_subspaces == 16;
    okay &= marked_serial.exact_distance_calls == 2;
    okay &= marked_parallel.supports_processed ==
            marked_serial.supports_processed;
    okay &= marked_parallel.eligible_supports ==
            marked_serial.eligible_supports;
    okay &= marked_parallel.marked_orbit_enumerated_supports ==
            marked_serial.marked_orbit_enumerated_supports;
    okay &= marked_parallel.nondegenerate_subspaces ==
            marked_serial.nondegenerate_subspaces;
    okay &= marked_parallel.marked_code_orbits ==
            marked_serial.marked_code_orbits;
    okay &= marked_parallel.exact_distance_calls ==
            marked_serial.exact_distance_calls;
    okay &= marked_parallel.raw_tensor_forms == marked_serial.raw_tensor_forms;
    okay &= marked_parallel.canonical_output_orbits ==
            marked_serial.canonical_output_orbits;
    okay &= marked_parallel.radical_dimension_counts ==
            marked_serial.radical_dimension_counts;
    okay &= marked_parallel.raw_tensor_keys == marked_serial.raw_tensor_keys;
    okay &= marked_parallel.canonical_output_keys ==
            marked_serial.canonical_output_keys;
    okay &= marked_parallel.positive_supports.size() ==
            marked_serial.positive_supports.size();
    for (std::size_t index = 0;
         index < std::min(
                     marked_parallel.positive_supports.size(),
                     marked_serial.positive_supports.size());
         ++index) {
      const auto& left = marked_parallel.positive_supports[index];
      const auto& right = marked_serial.positive_supports[index];
      okay &= left.support_offset == right.support_offset;
      okay &= left.record_index == right.record_index;
      okay &= left.output_keys == right.output_keys;
    }
    okay &= marked_parallel.pareto_witnesses.size() ==
            marked_serial.pareto_witnesses.size();
    for (std::size_t index = 0;
         index < std::min(
                     marked_parallel.pareto_witnesses.size(),
                     marked_serial.pareto_witnesses.size());
         ++index) {
      const auto& left = marked_parallel.pareto_witnesses[index];
      const auto& right = marked_serial.pareto_witnesses[index];
      okay &= left.output == right.output;
      okay &= left.known_class_id == right.known_class_id;
      okay &= left.distance == right.distance;
      okay &= left.error_coefficient == right.error_coefficient;
      okay &= left.protocol_length == right.protocol_length;
      okay &= left.space_footprint == right.space_footprint;
      okay &= left.record_index == right.record_index;
      okay &= left.source_index == right.source_index;
      okay &= left.source_id == right.source_id;
      okay &= left.source_family == right.source_family;
      okay &= left.source_space_length == right.source_space_length;
      okay &= left.ambient_dimension == right.ambient_dimension;
      okay &= left.parity_case == right.parity_case;
      okay &= left.origin == right.origin;
      okay &= left.points == right.points;
      okay &= left.canonical_logical_rows == right.canonical_logical_rows;
    }

  }

  std::vector<std::uint32_t> cube(16);
  for (std::uint32_t point = 0; point < cube.size(); ++point) cube[point] = point;
  const std::array<utsp::Mask, 1> labels{(utsp::Mask{1} << 16) - 2};
  const auto distance =
      utsp::z_distance_and_error_coefficient(cube, 4, labels);
  okay &= distance.distance == 3 && distance.error_coefficient == 35;
  const auto cube_distance_four =
      utsp::build_logical_label_space(cube, 4, true, 4);
  okay &= cube_distance_four.quotient_dimension() == 0;

  const std::array<utsp::ProtocolSupport, 1> catalogue_supports{
      utsp::ProtocolSupport{
          0,
          0,
          "self_test_cube",
          "self_test",
          16,
          4,
          0,
          0,
          cube,
      }};
  const auto catalogue = utsp::classify_protocol_frontier(
      catalogue_supports,
      1,
      3,
      24,
      utsp::ProtocolEnumerationMode::raw_subspaces,
      14,
      1,
      1,
      0,
      0,
      std::numeric_limits<std::uint32_t>::max(),
      true);
  okay &= catalogue.canonical_output_orbits == 1;
  okay &= catalogue.pareto_witnesses.size() == 1;
  okay &= catalogue.positive_supports.size() == 1;
  if (catalogue.positive_supports.size() == 1) {
    const auto& positive = catalogue.positive_supports.front();
    okay &= positive.support_offset == 0 && positive.record_index == 0;
    okay &= positive.output_keys.size() == 1;
  }
  if (catalogue.pareto_witnesses.size() == 1) {
    const auto& witness = catalogue.pareto_witnesses.front();
    okay &= witness.known_class_id == "D1_01";
    okay &= witness.protocol_length == 16;
    okay &= witness.space_footprint == 5;
    okay &= witness.distance == 3;
    okay &= witness.error_coefficient == 35;
  }
  const std::array<utsp::Mask, 3> independent_t_rows{1, 2, 4};
  const std::array<utsp::Mask, 4> one_idle_row{1, 2, 4, 0};
  okay &= utsp::intrinsic_tensor_dimension(independent_t_rows) == 3;
  okay &= utsp::intrinsic_tensor_dimension(one_idle_row) == 3;
  const auto radical_matches_definition = [](
                                               std::uint32_t q,
                                               std::uint64_t signature) {
    const auto packed_basis =
        utsp::cubic_tensor_radical_basis_from_signature(q, signature);
    std::vector<utsp::Mask> radical_basis(
        packed_basis.begin(), packed_basis.end());
    const auto vector_count = std::uint32_t{1} << q;
    for (std::uint32_t vector = 0; vector < vector_count; ++vector) {
      bool direct_radical = true;
      for (std::uint32_t left = 0; left < q && direct_radical; ++left) {
        for (std::uint32_t right = 0; right < q; ++right) {
          if (utsp::cubic_tensor_value(
                  q,
                  signature,
                  vector,
                  std::uint32_t{1} << left,
                  std::uint32_t{1} << right)) {
            direct_radical = false;
            break;
          }
        }
      }
      if ((utsp::reduce_vector(vector, radical_basis) == 0) !=
          direct_radical) {
        return false;
      }
    }
    return true;
  };
  for (std::uint32_t q = 1; q <= 3; ++q) {
    const auto signature_count =
        std::uint64_t{1} << utsp::cubic_tensor_term_count(q);
    for (std::uint64_t signature = 0; signature < signature_count;
         ++signature) {
      okay &= radical_matches_definition(q, signature);
    }
  }
  for (const auto q : {4U, 5U}) {
    const auto mask =
        (std::uint64_t{1} << utsp::cubic_tensor_term_count(q)) - 1;
    for (std::uint64_t sample = 0; sample < 1024; ++sample) {
      okay &= radical_matches_definition(
          q, splitmix64(sample ^ (std::uint64_t{q} << 32U)) & mask);
    }
  }
  std::array<utsp::Mask, 3> ccz_rows{};
  for (std::uint32_t point = 1; point < 8; ++point) {
    for (std::uint32_t coordinate = 0; coordinate < 3; ++coordinate) {
      if ((point >> coordinate) & 1U) {
        ccz_rows[coordinate] |= utsp::Mask{1} << (point - 1);
      }
    }
  }
  okay &= utsp::tensor_bilinear_rank(ccz_rows) == 0;
  okay &= utsp::tensor_radical_basis(ccz_rows).empty();
  okay &= utsp::nondegenerate_tensor_hyperplane_count(ccz_rows) == 0;
  okay &= utsp::is_primitive_exceptional_tensor(ccz_rows);
  const auto ccz_signature = utsp::cubic_tensor_word(ccz_rows);
  okay &= ccz_signature == 0x40;
  okay &= utsp::cubic_tensor_radical_basis_from_signature(
              3, ccz_signature)
              .empty();
  okay &= utsp::nondegenerate_tensor_hyperplane_count_from_signature(
              3, ccz_signature) == 0;
  const auto q3_zero_primitive = utsp::census_primitive_tensor_sector(
      3, utsp::PrimitiveGramSector::zero);
  okay &= q3_zero_primitive.tensor_count == 2;
  okay &= q3_zero_primitive.tensor_orbit_count == 2;
  okay &= q3_zero_primitive.primitive_tensor_count == 1;
  okay &= q3_zero_primitive.primitive_tensor_orbit_count == 1;
  const std::array<utsp::Mask, 5> primitive_q5_rows{
      ccz_rows[0] | (ccz_rows[0] << 7U),
      ccz_rows[1],
      ccz_rows[2],
      ccz_rows[1] << 7U,
      ccz_rows[2] << 7U,
  };
  okay &= utsp::tensor_bilinear_rank(primitive_q5_rows) == 0;
  okay &= utsp::tensor_radical_basis(primitive_q5_rows).empty();
  okay &=
      utsp::nondegenerate_tensor_hyperplane_count(primitive_q5_rows) == 0;
  okay &= utsp::is_primitive_exceptional_tensor(primitive_q5_rows);
  const std::array<utsp::Mask, 5> rank_one_primitive_q5_rows{
      0x1555, 0x1806, 0x198, 0x1e0, 0x1e00};
  okay &= utsp::cubic_tensor_word(rank_one_primitive_q5_rows) == 0x60001;
  okay &=
      utsp::is_primitive_exceptional_tensor(rank_one_primitive_q5_rows);
  const std::array<utsp::Mask, 6> decomposable_q6_rows{
      ccz_rows[0],
      ccz_rows[1],
      ccz_rows[2],
      ccz_rows[0] << 7U,
      ccz_rows[1] << 7U,
      ccz_rows[2] << 7U,
  };
  okay &= utsp::tensor_radical_basis(decomposable_q6_rows).empty();
  okay &=
      utsp::nondegenerate_tensor_hyperplane_count(decomposable_q6_rows) == 49;
  okay &= !utsp::is_primitive_exceptional_tensor(decomposable_q6_rows);
  const auto q7_authority =
      utsp::census_q7_primitive_tensors_from_authority(1);
  okay &= q7_authority.alternating_orbits.size() == 12;
  okay &= q7_authority.alternating_orbit_mass == (std::uint64_t{1} << 35);
  okay &= q7_authority.rank_one_pairs_tested == 1524;
  okay &= q7_authority.rank_one_nondegenerate_pairs == 954;
  okay &= q7_authority.rank_one_primitive_pairs == 1;
  okay &= q7_authority.primitive_orbits.size() == 2;
  if (q7_authority.primitive_orbits.size() == 2) {
    okay &= q7_authority.primitive_orbits[0].standard_tensor_signature ==
            0x42010000000ULL;
    okay &= q7_authority.primitive_orbits[1].standard_tensor_signature ==
            0x42010000001ULL;
  }
  const auto q8_authority =
      utsp::screen_q8_primitive_tensors_from_authority();
  okay &= q8_authority.alternating_orbits.size() == 32;
  okay &= q8_authority.alternating_orbit_mass == (std::uint64_t{1} << 56);
  okay &= q8_authority.rank_one_pairs_tested == 8160;
  okay &= q8_authority.rank_one_nondegenerate_pairs == 5868;
  okay &= q8_authority.rank_one_primitive_pairs == 0;
  okay &= q8_authority.primitive_candidates.empty();
  okay &= q8_authority.full_alternating_sector_theorem_excluded;
  const auto hard_rank_one_sector =
      utsp::rank_at_most_one_tensor_sector(hard_label_space);
  const auto hard_even_sector =
      utsp::even_weight_tensor_sector(hard_label_space);
  okay &= hard_rank_one_sector.bilinear_form_rows.size() == 7;
  okay &= hard_even_sector.quotient_dimension() == 20;
  const std::array<std::uint32_t, 4> distance_four_points{1, 2, 4, 7};
  const std::array<utsp::Mask, 1> distance_four_label{1};
  const auto distance_four =
      utsp::DistanceContext(distance_four_points, 3).evaluate(distance_four_label);
  okay &= distance_four.distance == 4 && distance_four.error_coefficient == 1;

  const utsp::KnownOutputClassifier classifier;
  const auto* t = classifier.classify(1, 0b1);
  const auto* cs = classifier.classify(2, 0b100);
  const auto* ccz = classifier.classify(3, 0b1000000);
  okay &= t && t->class_id == "D1_01";
  okay &= cs && cs->class_id == "D2_02";
  okay &= ccz && ccz->class_id == "D3_06";
  const auto cs_completion =
      classifier.minimum_feasible_completion(2, 0b100, 12, false, 40);
  okay &= cs_completion &&
          cs_completion->columns == std::vector<std::uint32_t>({1, 2, 3});
  utsp::SparseCompletionDecoder cs_decoder(2);
  const auto sparse_cs = cs_decoder.decode(0b100);
  okay &= sparse_cs && sparse_cs->columns == std::vector<std::uint32_t>({1, 2, 3}) &&
          sparse_cs->canonical_form.points ==
              std::vector<std::uint32_t>({1, 2, 3});
  utsp::SparseCompletionDecoder ccz_decoder(3);
  const auto sparse_ccz = ccz_decoder.decode(0b1000000);
  okay &= sparse_ccz && sparse_ccz->weight() == 7;
  const auto canonical_cs = utsp::canonicalize_cubic_tensor_direct(2, 0b100);
  okay &= canonical_cs.canonical_key == utsp::interleaved_tensor_key(
                                             2,
                                             0b100,
                                             canonical_cs
                                                 .new_basis_in_old_coordinates);
  const std::array<std::uint64_t, 2> cs_rows{0b101, 0b110};
  const auto canonical_cs_labels =
      utsp::canonicalize_cubic_tensor_from_labels(cs_rows);
  okay &= canonical_cs_labels.canonical_key == canonical_cs.canonical_key;
  const auto t3_centroid = utsp::analyze_tensor_centroid(independent_t_rows);
  okay &= t3_centroid.centroid_basis.size() == 3 && t3_centroid.split.has_value();
  const auto t3_decomposition =
      utsp::canonicalize_tensor_by_decomposition(independent_t_rows);
  okay &= t3_decomposition.components.size() == 3 &&
          t3_decomposition.radical_basis.empty();
  const std::array<std::uint32_t, 6> puncture_support{0, 1, 2, 4, 8, 12};
  const std::array<std::uint32_t, 3> puncture_template{0, 1, 2};
  const auto punctures = utsp::find_hereditary_punctures(
      puncture_support, 4, puncture_template, 2);
  okay &= punctures.complete &&
          std::find(
              punctures.puncture_index_masks.begin(),
              punctures.puncture_index_masks.end(),
              std::uint64_t{0b111}) != punctures.puncture_index_masks.end();
  const std::array<std::uint32_t, 3> affine_triangle{0, 1, 2};
  const auto triangle_orbits = utsp::affine_origin_orbits(affine_triangle, 2);
  std::vector<std::size_t> triangle_orbit_sizes;
  for (const auto& orbit : triangle_orbits.orbits) {
    triangle_orbit_sizes.push_back(orbit.members.size());
  }
  std::sort(triangle_orbit_sizes.begin(), triangle_orbit_sizes.end());
  okay &= triangle_orbits.complete &&
          triangle_orbit_sizes == std::vector<std::size_t>({1, 3});

  std::cout << "{\"schema\":\"utsp-native-self-test-v1\",\"status\":\""
            << (okay ? "pass" : "fail") << "\"}\n";
  return okay ? 0 : 1;
}

int inspect_support_automorphisms(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  const auto label_space =
      utsp::build_logical_label_space(points, ambient, false);
  const auto subgroup = utsp::elementary_support_automorphisms(label_space);
  std::cout << "{\"schema\":\"utsp-support-automorphisms-v1\","
            << "\"ambient_dimension\":" << subgroup.ambient_dimension
            << ",\"quotient_dimension\":" << subgroup.quotient_dimension
            << ",\"elementary_generator_count\":"
            << subgroup.elementary_generator_count
            << ",\"subgroup_order\":" << subgroup.subgroup_order
            << ",\"generator_count\":" << subgroup.generators.size()
            << ",\"generators\":[";
  bool first_generator = true;
  for (const auto& generator : subgroup.generators) {
    if (!first_generator) std::cout << ',';
    first_generator = false;
    std::cout << "{\"ambient_basis_images\":[";
    for (std::size_t index = 0;
         index < generator.ambient_basis_images.size(); ++index) {
      if (index != 0) std::cout << ',';
      std::cout << generator.ambient_basis_images[index];
    }
    std::cout << "],\"quotient_basis_images\":[";
    for (std::size_t index = 0;
         index < generator.quotient_basis_images.size(); ++index) {
      if (index != 0) std::cout << ',';
      std::cout << '\"' << utsp::hex_mask(
          generator.quotient_basis_images[index]) << '\"';
    }
    std::cout << "]}";
  }
  std::cout << "]}\n";
  return 0;
}

int inspect_support_vector_orbits(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  const auto label_space =
      utsp::build_logical_label_space(points, ambient, false);
  const auto subgroup = utsp::elementary_support_automorphisms(label_space);
  const auto orbits = utsp::quotient_vector_orbits(subgroup);
  std::uint64_t covered = 1;
  std::cout << "{\"schema\":\"utsp-support-vector-orbits-v1\","
            << "\"ambient_dimension\":" << subgroup.ambient_dimension
            << ",\"quotient_dimension\":" << subgroup.quotient_dimension
            << ",\"elementary_generator_count\":"
            << subgroup.elementary_generator_count
            << ",\"subgroup_order\":" << subgroup.subgroup_order
            << ",\"generator_count\":" << subgroup.generators.size()
            << ",\"orbits\":[";
  for (std::size_t index = 0; index < orbits.size(); ++index) {
    if (index != 0) std::cout << ',';
    covered += orbits[index].size;
    std::cout << "{\"representative\":\""
              << utsp::hex_mask(orbits[index].representative)
              << "\",\"size\":" << orbits[index].size << '}';
  }
  std::cout << "],\"covered_vectors_including_zero\":" << covered
            << "}\n";
  return 0;
}

int inspect_isotropic_orbit_cover(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto depth_text = argument_value(argc, argv, "--orbit-depth", false);
  const auto orbit_depth = depth_text.empty()
                               ? std::uint32_t{2}
                               : static_cast<std::uint32_t>(
                                     parse_unsigned(depth_text));
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  const auto label_space =
      utsp::build_logical_label_space(points, ambient, false);
  const auto subgroup = utsp::elementary_support_automorphisms(label_space);
  std::map<std::uint32_t, std::uint64_t> radical_dimensions;
  std::uint64_t nondegenerate = 0;
  std::uint64_t checksum = 0;
  std::uint64_t candidate_index = 0;
  const auto summary = utsp::enumerate_isotropic_orbit_cover(
      label_space,
      q,
      subgroup,
      orbit_depth,
      [&](std::span<const utsp::Mask> quotient_rows,
          std::span<const utsp::Mask> label_rows) {
        (void)quotient_rows;
        const auto radical = static_cast<std::uint32_t>(
            utsp::tensor_radical_basis(label_rows).size());
        ++radical_dimensions[radical];
        if (radical == 0) ++nondegenerate;
        std::uint64_t tensor_hash = 0;
        for (const auto word : utsp::cubic_tensor_words(label_rows)) {
          tensor_hash = splitmix64(tensor_hash ^ word);
        }
        checksum ^= splitmix64(tensor_hash ^ splitmix64(candidate_index++));
      });
  std::cout << "{\"schema\":\"utsp-isotropic-orbit-cover-v1\","
            << "\"ambient_dimension\":" << ambient
            << ",\"quotient_dimension\":"
            << label_space.quotient_dimension()
            << ",\"logical_dimension\":" << q
            << ",\"orbit_depth\":" << orbit_depth
            << ",\"generator_count\":" << subgroup.generators.size()
            << ",\"vector_orbits\":" << summary.vector_orbits
            << ",\"candidate_subspaces\":"
            << summary.candidate_subspaces
            << ",\"nondegenerate_candidates\":" << nondegenerate
            << ",\"radical_dimensions\":{";
  bool first = true;
  for (const auto& [dimension, count] : radical_dimensions) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '\"' << dimension << "\":" << count;
  }
  std::cout << "},\"tensor_checksum\":\"" << std::hex << checksum
            << std::dec << "\"}\n";
  return 0;
}

int inspect_canonical_isotropic_orbits(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto target_dimension_text =
      argument_value(argc, argv, "--target-dimension", false);
  const auto target_dimension = target_dimension_text.empty()
                                    ? q
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(
                                              target_dimension_text));
  const auto target_seed_text =
      argument_value(argc, argv, "--target-nondegenerate-seed", false);
  const auto target_seed = target_seed_text.empty()
                               ? std::uint32_t{0}
                               : static_cast<std::uint32_t>(
                                     parse_unsigned(target_seed_text));
  const auto emit_text = argument_value(argc, argv, "--emit", false);
  const auto emit = emit_text.empty()
                        ? std::uint64_t{0}
                        : parse_unsigned(emit_text);
  const auto emit_orbit_data = has_flag(argc, argv, "--emit-orbit-data");
  if (emit_orbit_data && emit == 0) {
    throw std::invalid_argument("--emit-orbit-data requires positive --emit");
  }
  const auto sector = argument_value(argc, argv, "--sector", false);
  const auto target_signatures_text =
      argument_value(argc, argv, "--target-signatures", false);
  const auto target_signatures = target_signatures_text.empty()
                                     ? std::vector<utsp::Mask>{}
                                     : parse_masks(target_signatures_text);
  const auto hereditary_nondegenerate =
      has_flag(argc, argv, "--hereditary-nondegenerate");
  const auto complete_nondegenerate =
      has_flag(argc, argv, "--complete-nondegenerate");
  const auto successor_nondegenerate =
      has_flag(argc, argv, "--successor-nondegenerate");
  const auto primitive_nondegenerate =
      has_flag(argc, argv, "--primitive-nondegenerate");
  const auto q5_zero_hyperplane_primitive =
      has_flag(argc, argv, "--q5-zero-hyperplane-primitive");
  const auto q5_zero_hyperplane_complete =
      has_flag(argc, argv, "--q5-zero-hyperplane-complete");
  const auto q5_q3_chain_cover =
      has_flag(argc, argv, "--q5-q3-chain-cover");
  const auto q5_q3_chain_profile_presence =
      has_flag(argc, argv, "--q5-q3-chain-profile-presence");
  const auto q5_q4_hitting_set =
      has_flag(argc, argv, "--q5-q4-hitting-set");
  const auto retain_from_text =
      argument_value(argc, argv, "--retain-nondegenerate-from", false);
  const auto selected_modes =
      static_cast<unsigned int>(hereditary_nondegenerate) +
      static_cast<unsigned int>(complete_nondegenerate) +
      static_cast<unsigned int>(successor_nondegenerate) +
      static_cast<unsigned int>(primitive_nondegenerate) +
      static_cast<unsigned int>(q5_zero_hyperplane_primitive) +
      static_cast<unsigned int>(q5_zero_hyperplane_complete) +
      static_cast<unsigned int>(q5_q3_chain_cover) +
      static_cast<unsigned int>(q5_q3_chain_profile_presence) +
      static_cast<unsigned int>(q5_q4_hitting_set) +
      static_cast<unsigned int>(!retain_from_text.empty());
  if (selected_modes > 1) {
    throw std::invalid_argument(
        "choose one nondegenerate orbit-search mode");
  }
  if ((q5_zero_hyperplane_primitive || q5_zero_hyperplane_complete ||
       q5_q3_chain_cover || q5_q3_chain_profile_presence ||
       q5_q4_hitting_set) &&
      q != 5) {
    throw std::invalid_argument(
        "specialized q5 modes require --q 5");
  }
  if (!target_signatures.empty() &&
      (hereditary_nondegenerate || complete_nondegenerate ||
       successor_nondegenerate || primitive_nondegenerate ||
       q5_zero_hyperplane_primitive || q5_zero_hyperplane_complete ||
       q5_q3_chain_cover || q5_q3_chain_profile_presence ||
       q5_q4_hitting_set ||
       !retain_from_text.empty() ||
       (!sector.empty() && sector != "all"))) {
    throw std::invalid_argument(
        "target signatures cannot be combined with sector or radical filters");
  }
  if (target_signatures.empty() && target_seed != 0) {
    throw std::invalid_argument(
        "a target nondegenerate seed requires target signatures");
  }
  const auto retain_nondegenerate_from = retain_from_text.empty()
                                             ? std::uint32_t{0}
                                             : static_cast<std::uint32_t>(
                                                   parse_unsigned(
                                                       retain_from_text));
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto workers = workers_text.empty()
                           ? std::uint32_t{1}
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  auto label_space = utsp::build_logical_label_space(points, ambient, false);
  if (sector.empty() || sector == "all") {
    // The complete simultaneous-isotropy problem.
  } else if (sector == "rank-at-most-one") {
    label_space = utsp::rank_at_most_one_tensor_sector(label_space);
  } else if (sector == "even") {
    label_space = utsp::even_weight_tensor_sector(label_space);
  } else {
    throw std::invalid_argument("unknown tensor sector: " + sector);
  }
  if (q5_q3_chain_profile_presence) {
    const auto started = Clock::now();
    const auto profile =
        utsp::enumerate_q5_q3_chain_output_profile_presence(
            label_space, workers);
    const auto finished = Clock::now();
    const auto& traversal = profile.traversal;
    std::cout
        << "{\"schema\":\"utsp-q5-output-profile-presence-v1\""
        << ",\"ambient_dimension\":" << ambient
        << ",\"quotient_dimension\":" << label_space.quotient_dimension()
        << ",\"logical_dimension\":5"
        << ",\"workers\":" << workers
        << ",\"parent_orbits_processed\":" << traversal.parent_orbits
        << ",\"extension_vector_orbits\":"
        << traversal.extension_vector_orbits
        << ",\"filtered_degenerate_children\":"
        << traversal.filtered_degenerate_children
        << ",\"filtered_target_children\":"
        << traversal.filtered_target_children
        << ",\"canonical_parent_rejections\":"
        << traversal.canonical_parent_rejections
        << ",\"duplicate_children\":" << traversal.duplicate_children
        << ",\"canonical_search_nodes\":"
        << traversal.canonical_search_nodes
        << ",\"support_automorphism_group_order\":"
        << traversal.support_automorphism_group_order
        << ",\"output_class_count\":"
        << profile.canonical_output_keys.size()
        << ",\"canonical_output_keys\":[";
    for (std::size_t index = 0;
         index < profile.canonical_output_keys.size(); ++index) {
      if (index != 0) std::cout << ',';
      std::cout << '[';
      for (std::size_t word = 0;
           word < profile.canonical_output_keys[index].size(); ++word) {
        if (word != 0) std::cout << ',';
        std::cout << '\"'
                  << utsp::hex_mask(profile.canonical_output_keys[index][word])
                  << '\"';
      }
      std::cout << ']';
    }
    std::cout << "],\"seconds\":" << std::setprecision(12)
              << seconds(started, finished) << "}\n";
    return 0;
  }
  const auto started = Clock::now();
  const auto level = complete_nondegenerate
                         ? utsp::enumerate_complete_nondegenerate_subspace_orbits(
                               label_space, q, workers)
                     : successor_nondegenerate
                         ? utsp::enumerate_nondegenerate_successor_subspace_orbits(
                               label_space, q, workers)
                     : q5_zero_hyperplane_complete
                         ? utsp::enumerate_complete_q5_zero_hyperplane_subspace_orbits(
                               label_space, workers)
                     : q5_zero_hyperplane_primitive
                         ? utsp::enumerate_q5_zero_hyperplane_primitive_subspace_orbits(
                               label_space, workers)
                     : q5_q3_chain_cover
                         ? utsp::enumerate_q5_q3_chain_cover_subspace_orbits(
                               label_space, workers)
                     : q5_q4_hitting_set
                         ? utsp::enumerate_q5_q4_hitting_set_subspace_orbits(
                               label_space, workers)
                     : primitive_nondegenerate
                         ? utsp::enumerate_primitive_subspace_orbits(
                               label_space, q, workers)
                     : !target_signatures.empty()
                         ? utsp::enumerate_target_tensor_subspace_orbits(
                               label_space,
                               q,
                               target_dimension,
                               target_signatures,
                               target_seed,
                               workers)
                     : hereditary_nondegenerate
                         ? utsp::enumerate_hereditary_nondegenerate_subspace_orbits(
                               label_space, q, workers)
                     : retain_nondegenerate_from != 0
                         ? utsp::enumerate_dimension_filtered_nondegenerate_subspace_orbits(
                               label_space,
                               q,
                               retain_nondegenerate_from,
                               workers)
                         : utsp::enumerate_isotropic_subspace_orbits(
                               label_space, q, workers);
  const auto finished = Clock::now();
  std::map<std::uint32_t, std::uint64_t> radical_dimensions;
  std::map<std::uint32_t, std::uint64_t> bilinear_ranks;
  std::uint64_t checksum = 1469598103934665603ULL;
  std::uint64_t canonical_key_checksum = 1469598103934665603ULL;
  std::uint64_t orbit_size_checksum = 1469598103934665603ULL;
  if (level.representatives.size() != level.canonical_keys.size() ||
      level.representatives.size() != level.orbit_sizes.size()) {
    throw std::logic_error("isotropic orbit output vectors are misaligned");
  }
  for (const auto& key : level.canonical_keys) {
    for (const auto word : key) {
      canonical_key_checksum ^= word;
      canonical_key_checksum *= 1099511628211ULL;
    }
  }
  for (const auto orbit_size : level.orbit_sizes) {
    orbit_size_checksum ^= orbit_size;
    orbit_size_checksum *= 1099511628211ULL;
  }
  for (const auto& quotient_rows : level.representatives) {
    std::vector<utsp::Mask> label_rows;
    label_rows.reserve(quotient_rows.size());
    for (const auto row : quotient_rows) {
      label_rows.push_back(
          utsp::linear_combination(row, label_space.quotient_basis));
      checksum ^= row;
      checksum *= 1099511628211ULL;
    }
    ++radical_dimensions[static_cast<std::uint32_t>(
        utsp::tensor_radical_basis(label_rows).size())];
    ++bilinear_ranks[utsp::tensor_bilinear_rank(label_rows)];
  }
  std::cout << "{\"schema\":\"utsp-canonical-isotropic-orbits-v1\"," 
            << "\"ambient_dimension\":" << ambient
            << ",\"quotient_dimension\":"
            << label_space.quotient_dimension()
            << ",\"logical_dimension\":" << q
            << ",\"workers\":" << workers
            << ",\"sector\":\""
            << (sector.empty() ? "all" : sector) << "\""
            << ",\"target_signature_count\":"
            << target_signatures.size()
            << ",\"target_dimension\":" << target_dimension
            << ",\"target_nondegenerate_seed\":" << target_seed
            << ",\"parent_orbits_processed\":" << level.parent_orbits
            << ",\"extension_vector_orbits\":"
            << level.extension_vector_orbits
            << ",\"hereditary_nondegenerate\":"
            << (hereditary_nondegenerate ? "true" : "false")
            << ",\"complete_nondegenerate\":"
            << (complete_nondegenerate ? "true" : "false")
            << ",\"successor_nondegenerate\":"
            << (successor_nondegenerate ? "true" : "false")
            << ",\"primitive_nondegenerate\":"
            << (primitive_nondegenerate ? "true" : "false")
            << ",\"q5_zero_hyperplane_primitive\":"
            << (q5_zero_hyperplane_primitive ? "true" : "false")
            << ",\"q5_zero_hyperplane_complete\":"
            << (q5_zero_hyperplane_complete ? "true" : "false")
            << ",\"q5_q3_chain_cover\":"
            << (q5_q3_chain_cover ? "true" : "false")
            << ",\"q5_q4_hitting_set\":"
            << (q5_q4_hitting_set ? "true" : "false")
            << ",\"retain_nondegenerate_from\":"
            << retain_nondegenerate_from
            << ",\"filtered_degenerate_children\":"
            << level.filtered_degenerate_children
            << ",\"filtered_nonprimitive_children\":"
            << level.filtered_nonprimitive_children
            << ",\"filtered_target_children\":"
            << level.filtered_target_children
            << ",\"used_complete_primitive_target_recognizer\":"
            << (level.used_complete_primitive_target_recognizer
                    ? "true"
                    : "false")
            << ",\"canonical_parent_rejections\":"
            << level.canonical_parent_rejections
            << ",\"duplicate_children\":" << level.duplicate_children
            << ",\"canonical_search_nodes\":"
            << level.canonical_search_nodes
            << ",\"support_automorphism_group_order\":"
            << level.support_automorphism_group_order
            << ",\"weighted_subspace_count\":"
            << level.weighted_subspace_count
            << ",\"subspace_orbits\":" << level.representatives.size()
            << ",\"radical_dimensions\":{";
  bool first = true;
  for (const auto& [dimension, count] : radical_dimensions) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '\"' << dimension << "\":" << count;
  }
  std::cout << "},\"tensor_bilinear_ranks\":{";
  first = true;
  for (const auto& [rank, count] : bilinear_ranks) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '\"' << rank << "\":" << count;
  }
  std::cout << "},\"representative_checksum\":\"" << std::hex
            << checksum << "\",\"canonical_key_checksum\":\""
            << canonical_key_checksum << "\",\"orbit_size_checksum\":\""
            << orbit_size_checksum << std::dec << "\",\"seconds\":"
            << std::setprecision(12) << seconds(started, finished)
            << ",\"representatives\":[";
  const auto emitted = std::min<std::uint64_t>(
      emit, level.representatives.size());
  for (std::size_t index = 0; index < emitted; ++index) {
    if (index != 0) std::cout << ',';
    std::cout << '[';
    for (std::size_t row = 0;
         row < level.representatives[index].size(); ++row) {
      if (row != 0) std::cout << ',';
      std::cout << '\"'
                << utsp::hex_mask(level.representatives[index][row])
                << '\"';
    }
    std::cout << ']';
  }
  std::cout << ']';
  if (emit_orbit_data) {
    std::cout << ",\"canonical_keys\":[";
    for (std::size_t index = 0; index < emitted; ++index) {
      if (index != 0) std::cout << ',';
      std::cout << '[';
      for (std::size_t word = 0;
           word < level.canonical_keys[index].size(); ++word) {
        if (word != 0) std::cout << ',';
        std::cout << '\"' << utsp::hex_mask(level.canonical_keys[index][word])
                  << '\"';
      }
      std::cout << ']';
    }
    std::cout << "],\"orbit_sizes\":[";
    for (std::size_t index = 0; index < emitted; ++index) {
      if (index != 0) std::cout << ',';
      std::cout << level.orbit_sizes[index];
    }
    std::cout << ']';
  }
  std::cout << "}\n";
  return 0;
}

int tensor_primitive_coverage(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto workers = workers_text.empty()
                           ? std::uint32_t{1}
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  if (q < 2) {
    throw std::invalid_argument("primitive coverage requires q >= 2");
  }
  const auto hereditary_from_text =
      argument_value(argc, argv, "--hereditary-from", false);
  const auto hereditary_from = hereditary_from_text.empty()
                                    ? q - 1
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(
                                              hereditary_from_text));
  if (hereditary_from == 0 || hereditary_from >= q) {
    throw std::invalid_argument(
        "primitive coverage requires 1 <= --hereditary-from < q");
  }
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  const auto base = utsp::build_logical_label_space(points, ambient, false);

  const auto hereditary_started = Clock::now();
  const auto hereditary =
      utsp::enumerate_dimension_filtered_nondegenerate_subspace_orbits(
          base, q, hereditary_from, workers);
  const auto hereditary_finished = Clock::now();
  std::set<std::vector<utsp::Mask>> hereditary_keys(
      hereditary.canonical_keys.begin(), hereditary.canonical_keys.end());
  if (hereditary_keys.size() != hereditary.representatives.size()) {
    throw std::logic_error("hereditary canonical keys are not unique");
  }

  const auto rank_one_space = utsp::rank_at_most_one_tensor_sector(base);
  const auto rank_one_started = Clock::now();
  const auto rank_one = utsp::enumerate_primitive_subspace_orbits(
      rank_one_space, q, workers);
  const auto rank_one_finished = Clock::now();

  std::set<std::vector<utsp::Mask>> rank_one_keys;
  for (std::size_t index = 0; index < rank_one.representatives.size(); ++index) {
    std::vector<utsp::Mask> labels;
    for (const auto row : rank_one.representatives[index]) {
      labels.push_back(
          utsp::linear_combination(row, rank_one_space.quotient_basis));
    }
    if (utsp::tensor_bilinear_rank(labels) > 1 ||
        !utsp::tensor_radical_basis(labels).empty()) {
      throw std::logic_error("rank-one primitive sector emitted an invalid row space");
    }
    rank_one_keys.insert(rank_one.canonical_keys[index]);
  }

  std::set<std::vector<utsp::Mask>> symplectic_keys;
  utsp::IsotropicOrbitLevel even;
  double even_seconds = 0;
  if ((q & 1U) == 0) {
    const auto even_space = utsp::even_weight_tensor_sector(base);
    const auto even_started = Clock::now();
    even = utsp::enumerate_primitive_subspace_orbits(
        even_space, q, workers);
    even_seconds = seconds(even_started, Clock::now());
    for (std::size_t index = 0; index < even.representatives.size(); ++index) {
      std::vector<utsp::Mask> labels;
      for (const auto row : even.representatives[index]) {
        labels.push_back(
            utsp::linear_combination(row, even_space.quotient_basis));
      }
      if (utsp::tensor_bilinear_rank(labels) == q) {
        if (!utsp::tensor_bilinear_is_alternating(labels) ||
            !utsp::tensor_radical_basis(labels).empty()) {
          throw std::logic_error(
              "full-alternating primitive sector emitted an invalid row space");
        }
        symplectic_keys.insert(even.canonical_keys[index]);
      }
    }
  }

  std::uint64_t new_rank_one = 0;
  std::uint64_t new_symplectic = 0;
  auto complete_keys = hereditary_keys;
  for (const auto& key : rank_one_keys) {
    if (complete_keys.insert(key).second) ++new_rank_one;
  }
  for (const auto& key : symplectic_keys) {
    if (complete_keys.insert(key).second) ++new_symplectic;
  }
  std::cout
      << "{\"schema\":\"utsp-tensor-primitive-coverage-v1\"," 
      << "\"ambient_dimension\":" << ambient
      << ",\"logical_dimension\":" << q
      << ",\"workers\":" << workers
      << ",\"hereditary_from\":" << hereditary_from
      << ",\"hereditary_orbits\":" << hereditary_keys.size()
      << ",\"rank_at_most_one_candidates\":" << rank_one_keys.size()
      << ",\"new_rank_at_most_one_orbits\":" << new_rank_one
      << ",\"full_alternating_candidates\":" << symplectic_keys.size()
      << ",\"new_full_alternating_orbits\":" << new_symplectic
      << ",\"complete_nondegenerate_orbits\":" << complete_keys.size()
      << ",\"hereditary_seconds\":" << std::setprecision(12)
      << seconds(hereditary_started, hereditary_finished)
      << ",\"rank_at_most_one_seconds\":"
      << seconds(rank_one_started, rank_one_finished)
      << ",\"even_sector_seconds\":" << even_seconds
      << "}\n";
  return 0;
}

int primitive_tensor_census(int argc, char** argv) {
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto sector_text = argument_value(argc, argv, "--sector", false);
  const auto emit_text = argument_value(argc, argv, "--emit", false);
  const auto emit = emit_text.empty()
                        ? std::uint64_t{0}
                        : parse_unsigned(emit_text);
  const auto direct_keys = has_flag(argc, argv, "--direct-keys");
  if (direct_keys && q > 5) {
    throw std::invalid_argument(
        "--direct-keys is intentionally limited to q <= 5");
  }

  auto sectors = utsp::primitive_gram_sectors(q);
  if (!sector_text.empty() && sector_text != "all") {
    std::optional<utsp::PrimitiveGramSector> selected;
    for (const auto sector : sectors) {
      if (utsp::primitive_gram_sector_name(sector) == sector_text) {
        selected = sector;
        break;
      }
    }
    if (!selected) {
      throw std::invalid_argument("invalid primitive Gram sector");
    }
    sectors = {*selected};
  }

  const auto all_started = Clock::now();
  std::uint64_t total_primitive_orbits = 0;
  std::uint64_t total_primitive_tensors = 0;
  std::set<std::uint64_t> direct_key_set;
  std::cout << "{\"schema\":\"utsp-primitive-tensor-census-v1\"," 
            << "\"logical_dimension\":" << q
            << ",\"sectors\":[";
  for (std::size_t sector_index = 0; sector_index < sectors.size();
       ++sector_index) {
    if (sector_index != 0) std::cout << ',';
    const auto started = Clock::now();
    const auto census =
        utsp::census_primitive_tensor_sector(q, sectors[sector_index]);
    const auto finished = Clock::now();
    total_primitive_orbits += census.primitive_tensor_orbit_count;
    total_primitive_tensors += census.primitive_tensor_count;
    std::cout << "{\"sector\":\""
              << utsp::primitive_gram_sector_name(census.sector) << "\""
              << ",\"repeated_signature\":\"" << std::hex
              << census.repeated_signature << std::dec << "\""
              << ",\"alternating_dimension\":"
              << census.alternating_dimension
              << ",\"action_generator_count\":"
              << census.action_generator_count
              << ",\"gram_stabilizer_order\":"
              << census.gram_stabilizer_order
              << ",\"generated_symplectic_transvection_labels\":"
              << census.generated_symplectic_transvection_labels
              << ",\"tensor_count\":" << census.tensor_count
              << ",\"tensor_orbit_count\":"
              << census.tensor_orbit_count
              << ",\"nondegenerate_tensor_count\":"
              << census.nondegenerate_tensor_count
              << ",\"nondegenerate_tensor_orbit_count\":"
              << census.nondegenerate_tensor_orbit_count
              << ",\"primitive_tensor_count\":"
              << census.primitive_tensor_count
              << ",\"primitive_tensor_orbit_count\":"
              << census.primitive_tensor_orbit_count
              << ",\"nondegenerate_hyperplane_histogram\":[";
    for (std::size_t index = 0;
         index < census.nondegenerate_hyperplane_histogram.size(); ++index) {
      if (index != 0) std::cout << ',';
      const auto& bucket = census.nondegenerate_hyperplane_histogram[index];
      std::cout << "{\"nondegenerate_hyperplanes\":"
                << bucket.nondegenerate_hyperplanes
                << ",\"orbit_count\":" << bucket.orbit_count
                << ",\"tensor_count\":" << bucket.tensor_count << '}';
    }
    std::cout << "],\"primitive_orbits\":[";
    const auto emitted = std::min<std::uint64_t>(
        emit, census.primitive_orbits.size());
    for (std::size_t index = 0; index < emitted; ++index) {
      if (index != 0) std::cout << ',';
      const auto& orbit = census.primitive_orbits[index];
      std::cout << "{\"alternating_representative\":\"" << std::hex
                << orbit.alternating_representative
                << "\",\"tensor_signature\":\""
                << orbit.standard_tensor_signature << std::dec << "\""
                << ",\"orbit_size\":" << orbit.orbit_size
                << ",\"tensor_stabilizer_order\":"
                << orbit.stabilizer_order;
      if (direct_keys) {
        const auto canonical = utsp::canonicalize_cubic_tensor_direct(
            q, orbit.standard_tensor_signature);
        if (!direct_key_set.insert(canonical.canonical_key).second) {
          throw std::logic_error(
              "distinct primitive Gram-stabilizer orbits share a direct key");
        }
        std::cout << ",\"direct_canonical_key\":\"" << std::hex
                  << canonical.canonical_key << std::dec << "\"";
      }
      std::cout << '}';
    }
    std::cout << "],\"seconds\":" << std::setprecision(12)
              << seconds(started, finished) << '}';
  }
  std::cout << "],\"primitive_tensor_orbits\":"
            << total_primitive_orbits
            << ",\"primitive_tensor_count\":"
            << total_primitive_tensors
            << ",\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(all_started, Clock::now()) << "}\n";
  return 0;
}

int q7_primitive_tensor_authority(int argc, char** argv) {
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto workers = workers_text.empty()
                           ? std::max(1U, std::thread::hardware_concurrency())
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  const auto started = Clock::now();
  const auto census =
      utsp::census_q7_primitive_tensors_from_authority(workers);
  std::uint64_t zero_primitive_orbits = 0;
  std::uint64_t rank_one_primitive_orbits = 0;
  for (const auto& orbit : census.primitive_orbits) {
    if (orbit.sector == utsp::PrimitiveGramSector::zero) {
      ++zero_primitive_orbits;
    } else if (orbit.sector == utsp::PrimitiveGramSector::rank_one) {
      ++rank_one_primitive_orbits;
    }
  }
  std::cout
      << "{\"schema\":\"utsp-q7-primitive-tensor-authority-v1\"," 
      << "\"logical_dimension\":7"
      << ",\"workers\":" << workers
      << ",\"general_linear_group_order\":"
      << census.general_linear_group_order
      << ",\"alternating_tensor_count\":"
      << census.alternating_tensor_count
      << ",\"alternating_orbit_mass\":"
      << census.alternating_orbit_mass
      << ",\"authority_orbit_count\":"
      << census.alternating_orbits.size()
      << ",\"alternating_orbits\":[";
  for (std::size_t index = 0; index < census.alternating_orbits.size();
       ++index) {
    if (index != 0) std::cout << ',';
    const auto& orbit = census.alternating_orbits[index];
    std::cout << "{\"authority_id\":" << orbit.authority_id
              << ",\"alternating_signature\":\"" << std::hex
              << orbit.alternating_signature << std::dec << "\""
              << ",\"stabilizer_order\":" << orbit.stabilizer_order
              << ",\"orbit_size\":" << orbit.orbit_size
              << ",\"radical_dimension\":" << orbit.radical_dimension
              << ",\"nondegenerate_hyperplanes\":"
              << orbit.nondegenerate_hyperplanes << '}';
  }
  std::cout << "]"
            << ",\"rank_one_pairs_tested\":"
            << census.rank_one_pairs_tested
            << ",\"rank_one_nondegenerate_pairs\":"
            << census.rank_one_nondegenerate_pairs
            << ",\"rank_one_primitive_pairs\":"
            << census.rank_one_primitive_pairs
            << ",\"zero_gram_primitive_orbits\":"
            << zero_primitive_orbits
            << ",\"rank_one_primitive_orbits\":"
            << rank_one_primitive_orbits
            << ",\"primitive_orbits\":[";
  for (std::size_t index = 0; index < census.primitive_orbits.size();
       ++index) {
    if (index != 0) std::cout << ',';
    const auto& orbit = census.primitive_orbits[index];
    std::cout << "{\"sector\":\""
              << utsp::primitive_gram_sector_name(orbit.sector) << "\""
              << ",\"authority_id\":" << orbit.authority_id
              << ",\"rank_one_functional\":"
              << orbit.rank_one_functional
              << ",\"tensor_signature\":\"" << std::hex
              << orbit.standard_tensor_signature << std::dec << "\""
              << ",\"direct_canonical_key_words\":[";
    for (std::size_t word = 0; word < orbit.direct_canonical_key.size();
         ++word) {
      if (word != 0) std::cout << ',';
      std::cout << '"' << std::hex << orbit.direct_canonical_key[word]
                << std::dec << '"';
    }
    std::cout << "]"
              << ",\"represented_authority_pairs\":"
              << orbit.represented_authority_pairs << '}';
  }
  std::cout << "]"
            << ",\"primitive_tensor_orbits\":"
            << census.primitive_orbits.size()
            << ",\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return 0;
}

int q8_primitive_tensor_authority_screen() {
  const auto started = Clock::now();
  const auto screen = utsp::screen_q8_primitive_tensors_from_authority();
  std::uint64_t zero_primitive_orbits = 0;
  std::uint64_t rank_one_primitive_candidates = 0;
  for (const auto& candidate : screen.primitive_candidates) {
    if (candidate.sector == utsp::PrimitiveGramSector::zero) {
      ++zero_primitive_orbits;
    } else if (candidate.sector == utsp::PrimitiveGramSector::rank_one) {
      ++rank_one_primitive_candidates;
    }
  }
  std::cout
      << "{\"schema\":\"utsp-q8-primitive-tensor-authority-screen-v1\"," 
      << "\"logical_dimension\":8"
      << ",\"general_linear_group_order\":"
      << screen.general_linear_group_order
      << ",\"alternating_tensor_count\":"
      << screen.alternating_tensor_count
      << ",\"alternating_orbit_mass\":"
      << screen.alternating_orbit_mass
      << ",\"authority_orbit_count\":"
      << screen.alternating_orbits.size()
      << ",\"full_alternating_sector_theorem_excluded\":"
      << (screen.full_alternating_sector_theorem_excluded ? "true" : "false")
      << ",\"alternating_orbits\":[";
  for (std::size_t index = 0; index < screen.alternating_orbits.size();
       ++index) {
    if (index != 0) std::cout << ',';
    const auto& orbit = screen.alternating_orbits[index];
    std::cout << "{\"authority_id\":" << orbit.authority_id
              << ",\"alternating_state\":\"" << std::hex
              << orbit.alternating_state << std::dec << "\""
              << ",\"stabilizer_order\":" << orbit.stabilizer_order
              << ",\"orbit_size\":" << orbit.orbit_size
              << ",\"radical_dimension\":"
              << orbit.radical_dimension << '}';
  }
  std::cout << "]"
            << ",\"rank_one_pairs_tested\":"
            << screen.rank_one_pairs_tested
            << ",\"rank_one_nondegenerate_pairs\":"
            << screen.rank_one_nondegenerate_pairs
            << ",\"rank_one_primitive_pairs\":"
            << screen.rank_one_primitive_pairs
            << ",\"zero_gram_primitive_orbits\":"
            << zero_primitive_orbits
            << ",\"rank_one_primitive_candidates\":"
            << rank_one_primitive_candidates
            << ",\"primitive_candidates\":[";
  for (std::size_t index = 0; index < screen.primitive_candidates.size();
       ++index) {
    if (index != 0) std::cout << ',';
    const auto& candidate = screen.primitive_candidates[index];
    std::cout << "{\"sector\":\""
              << utsp::primitive_gram_sector_name(candidate.sector) << "\""
              << ",\"authority_id\":" << candidate.authority_id
              << ",\"rank_one_functional\":"
              << candidate.rank_one_functional
              << ",\"alternating_state\":\"" << std::hex
              << candidate.alternating_state << std::dec << "\"}";
  }
  std::cout << "]"
            << ",\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return 0;
}

int inspect_support(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto distance_text = argument_value(argc, argv, "--minimum-distance", false);
  const auto minimum_distance = distance_text.empty()
                                    ? std::uint32_t{3}
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(distance_text));
  const auto emit_text = argument_value(argc, argv, "--emit", false);
  const auto emit = emit_text.empty() ? std::uint64_t{0} : parse_unsigned(emit_text);
  const auto points = parse_points(argument_value(argc, argv, "--points"));
  const auto label_space =
      utsp::build_logical_label_space(points, ambient, true, minimum_distance);

  std::uint64_t emitted = 0;
  std::ostringstream records;
  std::uint64_t count = 0;
  if (q != 0) {
    count = utsp::enumerate_totally_isotropic_subspaces(
        label_space, q,
        [&](std::span<const utsp::Mask> quotient, std::span<const utsp::Mask> labels) {
        if (emitted >= emit) return;
        if (emitted) records << ',';
        records << "{\"quotient_rows\":[";
        for (std::size_t index = 0; index < quotient.size(); ++index) {
          if (index) records << ',';
          records << '"' << utsp::hex_mask(quotient[index]) << '"';
        }
        records << "],\"label_rows\":[";
        for (std::size_t index = 0; index < labels.size(); ++index) {
          if (index) records << ',';
          records << '"' << utsp::hex_mask(labels[index]) << '"';
        }
        records << "],\"tensor_words\":[";
        const auto tensor = utsp::cubic_tensor_words(labels);
        for (std::size_t index = 0; index < tensor.size(); ++index) {
          if (index) records << ',';
          records << '"' << utsp::hex_mask(tensor[index]) << '"';
        }
        records << "]}";
          ++emitted;
        });
  }

  std::cout << "{\"schema\":\"utsp-native-inspect-v1\",\"ambient_dimension\":"
            << ambient << ",\"support_length\":" << points.size()
            << ",\"minimum_distance\":" << minimum_distance
            << ",\"quotient_dimension\":" << label_space.quotient_dimension()
            << ",\"coordinate_masks\":";
  print_mask_array(label_space.coordinate_masks);
  std::cout << ",\"bilinear_form_rows\":[";
  for (std::size_t index = 0; index < label_space.bilinear_form_rows.size(); ++index) {
    if (index) std::cout << ',';
    print_mask_array(label_space.bilinear_form_rows[index]);
  }
  std::cout << ']'
            << ",\"stabilizer_basis\":";
  print_mask_array(label_space.stabilizer_basis);
  std::cout << ",\"quotient_basis\":";
  print_mask_array(label_space.quotient_basis);
  std::cout << ",\"isotropic_subspace_count\":" << count
            << ",\"subspaces\":[" << records.str() << "]}\n";
  return 0;
}

struct SourceSchurQuotientCache {
  bool valid = false;
  std::uint32_t dimension = 0;
  std::vector<std::uint32_t> points;
  std::vector<utsp::Mask> complement_basis;
  std::vector<utsp::Mask> pointed_basis;
  std::array<std::uint64_t, 64> feature_rows{};
  std::array<utsp::Mask, 64> feature_coefficients{};
};

[[nodiscard]] utsp::Mask support_width_mask(std::size_t width) {
  if (width > 64) throw std::invalid_argument("support width exceeds 64");
  return width == 64 ? ~utsp::Mask{0}
                     : (utsp::Mask{1} << width) - utsp::Mask{1};
}

[[nodiscard]] bool mask_dot(utsp::Mask left, utsp::Mask right) {
  return (std::popcount(left & right) & 1U) != 0;
}

[[nodiscard]] std::uint64_t quadratic_feature_vector(
    std::uint32_t point,
    std::uint32_t dimension) {
  const auto feature_count = 1U + dimension + dimension * (dimension - 1U) / 2U;
  if (feature_count > 64) {
    throw std::invalid_argument("quadratic feature vector exceeds uint64_t");
  }
  std::uint64_t result = 1;
  std::uint32_t feature = 1;
  for (std::uint32_t coordinate = 0; coordinate < dimension;
       ++coordinate, ++feature) {
    if ((point >> coordinate) & 1U) result |= std::uint64_t{1} << feature;
  }
  for (std::uint32_t left = 0; left < dimension; ++left) {
    for (std::uint32_t right = left + 1; right < dimension;
         ++right, ++feature) {
      if (((point >> left) & 1U) && ((point >> right) & 1U)) {
        result |= std::uint64_t{1} << feature;
      }
    }
  }
  return result;
}

[[nodiscard]] std::vector<std::uint32_t> recover_source_points(
    const PointedRecord& record) {
  const auto dimension = record.base_affine_dimension;
  if (dimension >= 31) {
    throw std::invalid_argument("source dimension exceeds uint32_t point format");
  }
  const auto low_mask = (std::uint32_t{1} << dimension) - 1U;
  std::vector<std::uint32_t> result;
  result.reserve(record.base_length);
  if (record.parity_case <= 2) {
    if (record.origin < 0) {
      throw std::invalid_argument("translated record has no source origin");
    }
    const auto origin = static_cast<std::uint32_t>(record.origin);
    for (const auto point : record.points) {
      if (record.parity_case == 2 && point == 0) continue;
      result.push_back(point ^ origin);
    }
    if (record.parity_case == 1) result.push_back(origin);
  } else if (record.parity_case == 3 || record.parity_case == 4) {
    const auto extension_bit = std::uint32_t{1} << dimension;
    for (const auto point : record.points) {
      if (record.parity_case == 4 && point == 0) continue;
      if ((point & extension_bit) == 0 || (point & ~low_mask) != extension_bit) {
        throw std::invalid_argument("hyperplane record has invalid source image");
      }
      result.push_back(point & low_mask);
    }
  } else {
    throw std::invalid_argument("unknown pointed-support parity case");
  }
  std::sort(result.begin(), result.end());
  if (result.size() != record.base_length ||
      std::adjacent_find(result.begin(), result.end()) != result.end()) {
    throw std::invalid_argument("failed to reconstruct source support");
  }
  return result;
}

[[nodiscard]] SourceSchurQuotientCache build_source_schur_cache(
    std::span<const std::uint32_t> points,
    std::uint32_t dimension) {
  if (points.empty() || points.size() > 64 || dimension >= 11) {
    throw std::invalid_argument("unsupported source for Schur quotient cache");
  }
  SourceSchurQuotientCache result;
  result.valid = true;
  result.dimension = dimension;
  result.points.assign(points.begin(), points.end());
  const auto one = support_width_mask(points.size());
  std::vector<utsp::Mask> coordinates(dimension, 0);
  for (std::size_t index = 0; index < points.size(); ++index) {
    if (points[index] >= (std::uint32_t{1} << dimension)) {
      throw std::invalid_argument("source point lies outside ambient space");
    }
    for (std::uint32_t coordinate = 0; coordinate < dimension; ++coordinate) {
      if ((points[index] >> coordinate) & 1U) {
        coordinates[coordinate] |= utsp::Mask{1} << index;
      }
    }
  }
  if (utsp::gf2_rank(coordinates) != dimension) {
    throw std::invalid_argument("source support does not span its ambient space");
  }

  std::vector<utsp::Mask> schur_rows{one};
  schur_rows.insert(schur_rows.end(), coordinates.begin(), coordinates.end());
  for (std::uint32_t left = 0; left < dimension; ++left) {
    for (std::uint32_t right = left + 1; right < dimension; ++right) {
      schur_rows.push_back(coordinates[left] & coordinates[right]);
    }
  }
  schur_rows = utsp::rref_basis(schur_rows);
  const auto dual = utsp::nullspace_basis(schur_rows, points.size());
  std::vector<utsp::Mask> affine_rows{one};
  affine_rows.insert(affine_rows.end(), coordinates.begin(), coordinates.end());
  auto span = utsp::rref_basis(affine_rows);
  if (span.size() != dimension + 1U) {
    throw std::invalid_argument("source affine rows are dependent");
  }
  for (const auto constraint : schur_rows) {
    for (const auto row : span) {
      if (mask_dot(constraint, row)) {
        throw std::invalid_argument("source is not unital triorthogonal");
      }
    }
  }
  for (const auto row : dual) {
    const auto residual = utsp::reduce_vector(row, span);
    if (residual == 0) continue;
    result.complement_basis.push_back(residual);
    span.push_back(residual);
    span = utsp::rref_basis(span);
  }
  result.pointed_basis = result.complement_basis;
  result.pointed_basis.push_back(one);

  for (std::size_t index = 0; index < points.size(); ++index) {
    auto value = quadratic_feature_vector(points[index], dimension);
    auto coefficients = utsp::Mask{1} << index;
    while (value != 0) {
      const auto pivot = 63U - static_cast<std::uint32_t>(std::countl_zero(value));
      if (result.feature_rows[pivot] == 0) {
        result.feature_rows[pivot] = value;
        result.feature_coefficients[pivot] = coefficients;
        break;
      }
      value ^= result.feature_rows[pivot];
      coefficients ^= result.feature_coefficients[pivot];
    }
  }
  return result;
}

[[nodiscard]] std::optional<utsp::Mask> source_feature_witness(
    const SourceSchurQuotientCache& cache,
    std::uint32_t point) {
  auto value = quadratic_feature_vector(point, cache.dimension);
  utsp::Mask coefficients = 0;
  while (value != 0) {
    const auto pivot = 63U - static_cast<std::uint32_t>(std::countl_zero(value));
    if (cache.feature_rows[pivot] == 0) return std::nullopt;
    value ^= cache.feature_rows[pivot];
    coefficients ^= cache.feature_coefficients[pivot];
  }
  return coefficients;
}

[[nodiscard]] utsp::Mask map_source_row(
    utsp::Mask row,
    const SourceSchurQuotientCache& cache,
    const PointedRecord& record,
    std::optional<std::size_t> skipped_source_index = std::nullopt) {
  utsp::Mask result = 0;
  for (std::size_t index = 0; index < cache.points.size(); ++index) {
    if (skipped_source_index == index || ((row >> index) & 1U) == 0) continue;
    const auto target = record.parity_case <= 2
                            ? cache.points[index] ^
                                  static_cast<std::uint32_t>(record.origin)
                            : cache.points[index] |
                                  (std::uint32_t{1} << cache.dimension);
    const auto found = std::lower_bound(
        record.points.begin(), record.points.end(), target);
    if (found == record.points.end() || *found != target) {
      throw std::invalid_argument("pointed support does not match its source map");
    }
    result |= utsp::Mask{1} << std::distance(record.points.begin(), found);
  }
  return result;
}

[[nodiscard]] std::vector<utsp::Mask> source_schur_pointed_quotient_basis(
    const SourceSchurQuotientCache& cache,
    const PointedRecord& record) {
  if (!cache.valid || record.base_affine_dimension != cache.dimension ||
      record.base_length != cache.points.size()) {
    throw std::invalid_argument("pointed record and source Schur cache disagree");
  }
  std::vector<utsp::Mask> source_basis;
  std::optional<std::size_t> skipped;
  if (record.parity_case <= 2) {
    if (record.origin < 0) {
      throw std::invalid_argument("translated pointed support has no origin");
    }
    const auto origin = static_cast<std::uint32_t>(record.origin);
    const auto found =
        std::lower_bound(cache.points.begin(), cache.points.end(), origin);
    const bool inside = found != cache.points.end() && *found == origin;
    if ((record.parity_case == 1 && !inside) ||
        (record.parity_case == 2 && inside)) {
      throw std::invalid_argument("pointed parity case disagrees with its origin");
    }
    source_basis = cache.pointed_basis;
    if (record.parity_case == 0 && inside) {
      const auto source_index =
          static_cast<std::size_t>(std::distance(cache.points.begin(), found));
      const auto origin_bit = utsp::Mask{1} << source_index;
      for (auto& row : source_basis) {
        if (row & origin_bit) row ^= origin_bit;
      }
    } else if (record.parity_case == 0 || record.parity_case == 2) {
      const auto witness = source_feature_witness(cache, origin);
      if (witness) source_basis.push_back(*witness);
    } else {
      skipped = static_cast<std::size_t>(
          std::distance(cache.points.begin(), found));
    }
  } else if (record.parity_case == 3 || record.parity_case == 4) {
    source_basis = cache.complement_basis;
  } else {
    throw std::invalid_argument("unknown pointed-support parity case");
  }

  std::vector<utsp::Mask> result;
  result.reserve(source_basis.size());
  for (const auto row : source_basis) {
    result.push_back(map_source_row(row, cache, record, skipped));
  }
  return result;
}

[[nodiscard]] std::vector<SourceSchurQuotientCache> build_source_schur_caches(
    const Manifest& manifest,
    std::uint32_t requested_workers) {
  std::vector<const PointedRecord*> exemplars(manifest.sources.size(), nullptr);
  for (const auto& record : manifest.supports) {
    if (record.source_index >= exemplars.size()) {
      throw std::invalid_argument("pointed record source index is out of range");
    }
    if (exemplars[record.source_index] == nullptr) {
      exemplars[record.source_index] = &record;
    }
  }
  std::vector<SourceSchurQuotientCache> result(manifest.sources.size());
  std::atomic<std::size_t> next{0};
  std::exception_ptr worker_error;
  std::mutex error_mutex;
  const auto worker_count =
      std::min<std::size_t>(requested_workers, manifest.sources.size());
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    workers.emplace_back([&] {
      try {
        while (true) {
          const auto source_index = next.fetch_add(1);
          if (source_index >= manifest.sources.size()) return;
          const auto* exemplar = exemplars[source_index];
          if (exemplar == nullptr) continue;
          const auto points = recover_source_points(*exemplar);
          result[source_index] = build_source_schur_cache(
              points, exemplar->base_affine_dimension);
        }
      } catch (...) {
        const std::lock_guard lock(error_mutex);
        if (!worker_error) worker_error = std::current_exception();
        next.store(manifest.sources.size());
      }
    });
  }
  for (auto& worker : workers) worker.join();
  if (worker_error) std::rethrow_exception(worker_error);
  return result;
}

int cache_manifest_label_spaces(int argc, char** argv) {
  namespace fs = std::filesystem;
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto positive_manifest_path =
      argument_value(argc, argv, "--positive-manifest", false);
  const bool positive_only = !positive_manifest_path.empty();
  const auto selected_dimension_text = argument_value(
      argc, argv, "--minimum-selected-quotient-dimension", false);
  const auto minimum_selected_quotient_dimension =
      selected_dimension_text.empty()
          ? std::uint32_t{1}
          : static_cast<std::uint32_t>(
                parse_unsigned(selected_dimension_text));
  const auto distance_text =
      argument_value(argc, argv, "--minimum-distance", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const bool source_schur = has_flag(argc, argv, "--source-schur");
  const auto minimum_distance = distance_text.empty()
                                    ? std::uint32_t{4}
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(distance_text));
  const auto requested_workers = workers_text.empty()
                                     ? std::uint32_t{1}
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  if (requested_workers == 0) {
    throw std::invalid_argument("label cache worker count must be positive");
  }
  if (minimum_distance < 3 || minimum_distance > 5) {
    throw std::invalid_argument("label cache supports minimum distance 3,4,5");
  }
  if (!positive_only && !selected_dimension_text.empty()) {
    throw std::invalid_argument(
        "selected quotient dimension requires --positive-manifest");
  }
  if (minimum_selected_quotient_dimension > 64) {
    throw std::invalid_argument(
        "selected quotient dimension exceeds 64");
  }

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "label cache requires an origin-expanded or origin-reduced manifest");
  }
  const auto read_finished = Clock::now();

  std::vector<SourceSchurQuotientCache> source_caches;
  double source_cache_seconds = 0;
  if (source_schur) {
    const auto source_cache_started = Clock::now();
    source_caches = build_source_schur_caches(manifest, requested_workers);
    const auto source_cache_finished = Clock::now();
    source_cache_seconds = seconds(
        source_cache_started, source_cache_finished);
  }

  const fs::path target(output_path);
  if (!target.parent_path().empty()) fs::create_directories(target.parent_path());
  const fs::path temporary = target.string() + ".tmp";
  std::map<std::uint32_t, std::uint64_t> dimension_histogram;
  Manifest positive_manifest;
  if (positive_only) {
    positive_manifest.format_version = 2;
    positive_manifest.length_limit = manifest.length_limit;
    positive_manifest.protocol_length_filter =
        manifest.protocol_length_filter;
    positive_manifest.base_support_only = false;
    positive_manifest.sources = manifest.sources;
  }
  const auto build_started = Clock::now();
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "could not open temporary label cache: " + temporary.string());
    }
    output.write(kLabelCacheMagic.data(), kLabelCacheMagic.size());
    write_little(output, std::uint32_t{1});
    write_little(output, minimum_distance);
    write_little(
        output,
        positive_only ? std::uint64_t{0}
                      : static_cast<std::uint64_t>(manifest.supports.size()));
    write_little(
        output,
        positive_only ? std::uint64_t{0}
                      : manifest_cache_fingerprint(manifest));

    constexpr std::size_t minimum_chunk = 4096;
    const auto chunk_size = std::max<std::size_t>(
        minimum_chunk, static_cast<std::size_t>(requested_workers) * 1024U);
    for (std::size_t chunk_begin = 0;
         chunk_begin < manifest.supports.size();
         chunk_begin += chunk_size) {
      const auto chunk_end =
          std::min(manifest.supports.size(), chunk_begin + chunk_size);
      const auto count = chunk_end - chunk_begin;
      std::vector<std::vector<utsp::Mask>> quotient_bases(count);
      std::atomic<std::size_t> next{0};
      const auto worker_count = std::min<std::size_t>(requested_workers, count);
      std::vector<std::exception_ptr> errors(worker_count);
      std::vector<std::thread> threads;
      threads.reserve(worker_count);
      for (std::size_t worker = 0; worker < worker_count; ++worker) {
        threads.emplace_back([&, worker] {
          while (true) {
            const auto offset = next.fetch_add(1);
            if (offset >= count) return;
            try {
              const auto& record = manifest.supports[chunk_begin + offset];
              if (source_schur) {
                if (record.source_index >= source_caches.size()) {
                  throw std::invalid_argument(
                      "pointed record source index is out of range");
                }
                auto basis = source_schur_pointed_quotient_basis(
                    source_caches[record.source_index], record);
                if (minimum_distance > 3) {
                  basis =
                      utsp::restrict_logical_label_quotient_basis_minimum_distance(
                          record.points,
                          record.ambient_dimension,
                          basis,
                          minimum_distance,
                          false);
                }
                quotient_bases[offset] = std::move(basis);
              } else {
                quotient_bases[offset] =
                    utsp::build_logical_label_quotient_basis(
                        record.points,
                        record.ambient_dimension,
                        false,
                        minimum_distance);
              }
            } catch (...) {
              errors[worker] = std::current_exception();
              next.store(count);
              return;
            }
          }
        });
      }
      for (auto& thread : threads) thread.join();
      for (const auto& error : errors) {
        if (error) std::rethrow_exception(error);
      }
      for (std::size_t offset = 0; offset < count; ++offset) {
        const auto& record = manifest.supports[chunk_begin + offset];
        const auto& basis = quotient_bases[offset];
        if (basis.size() > 64) {
          throw std::logic_error("label quotient dimension exceeds 64");
        }
        ++dimension_histogram[static_cast<std::uint32_t>(basis.size())];
        if (positive_only &&
            basis.size() < minimum_selected_quotient_dimension) {
          continue;
        }
        if (positive_only) positive_manifest.supports.push_back(record);
        write_little(output, static_cast<std::uint8_t>(basis.size()));
        for (const auto row : basis) write_little(output, row);
      }
    }
    if (positive_only) {
      output.seekp(16, std::ios::beg);
      write_little(
          output,
          static_cast<std::uint64_t>(positive_manifest.supports.size()));
      write_little(output, manifest_cache_fingerprint(positive_manifest));
      output.seekp(0, std::ios::end);
    }
    output.flush();
    if (!output) throw std::runtime_error("failed to flush temporary label cache");
  }
  const auto build_finished = Clock::now();
  if (positive_only) {
    write_manifest(positive_manifest, positive_manifest_path);
  }
  std::error_code error;
  fs::rename(temporary, target, error);
  if (error) {
    fs::remove(temporary);
    throw std::runtime_error("could not promote label cache: " + error.message());
  }

  std::cout << "{\"schema\":\"utsp-native-label-cache-summary-v1\"," 
            << "\"status\":\"complete\"," 
            << "\"minimum_distance\":" << minimum_distance << ','
            << "\"input_supports\":" << manifest.supports.size() << ','
            << "\"supports\":"
            << (positive_only ? positive_manifest.supports.size()
                              : manifest.supports.size())
            << ','
            << "\"positive_only\":"
            << (positive_only ? "true" : "false") << ','
            << "\"minimum_selected_quotient_dimension\":"
            << (positive_only ? minimum_selected_quotient_dimension : 0U)
            << ','
             << "\"workers\":" << requested_workers << ','
            << "\"source_schur\":"
            << (source_schur ? "true" : "false") << ','
            << "\"source_cache_seconds\":" << std::setprecision(12)
            << source_cache_seconds << ','
            << "\"cache_bytes\":" << fs::file_size(target) << ','
            << "\"positive_manifest_bytes\":";
  if (positive_only) {
    std::cout << fs::file_size(positive_manifest_path);
  } else {
    std::cout << "null";
  }
  std::cout << ','
            << "\"quotient_dimension_histogram\":[";
  bool first = true;
  for (const auto& [dimension, count] : dimension_histogram) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << dimension << ',' << count << ']';
  }
  std::cout << "],\"read_seconds\":" << std::setprecision(12)
            << seconds(read_started, read_finished)
            << ",\"build_seconds\":"
            << seconds(build_started, build_finished) << "}\n";
  return 0;
}

int cache_manifest_initial_label_spaces(int argc, char** argv) {
  namespace fs = std::filesystem;
  const auto input_path = argument_value(argc, argv, "--input");
  const auto distance_three_cache_path =
      argument_value(argc, argv, "--distance-three-cache");
  const auto distance_three_manifest_path =
      argument_value(argc, argv, "--distance-three-manifest");
  const auto distance_four_cache_path =
      argument_value(argc, argv, "--distance-four-cache");
  const auto distance_four_manifest_path =
      argument_value(argc, argv, "--distance-four-manifest");
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const bool source_schur = has_flag(argc, argv, "--source-schur");
  const auto requested_workers = workers_text.empty()
                                     ? std::uint32_t{1}
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  if (requested_workers == 0) {
    throw std::invalid_argument("initial label cache worker count must be positive");
  }
  const std::set<std::string> output_paths{
      distance_three_cache_path,
      distance_three_manifest_path,
      distance_four_cache_path,
      distance_four_manifest_path,
  };
  if (output_paths.size() != 4) {
    throw std::invalid_argument("initial label cache output paths must be distinct");
  }

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "initial label cache requires an origin-expanded or origin-reduced manifest");
  }
  const auto read_finished = Clock::now();

  std::vector<SourceSchurQuotientCache> source_caches;
  double source_cache_seconds = 0;
  if (source_schur) {
    const auto source_cache_started = Clock::now();
    source_caches = build_source_schur_caches(manifest, requested_workers);
    const auto source_cache_finished = Clock::now();
    source_cache_seconds = seconds(
        source_cache_started, source_cache_finished);
  }

  const fs::path distance_three_target(distance_three_cache_path);
  const fs::path distance_four_target(distance_four_cache_path);
  const fs::path distance_three_temporary =
      distance_three_target.string() + ".tmp";
  const fs::path distance_four_temporary =
      distance_four_target.string() + ".tmp";
  const std::array<fs::path, 4> initial_output_paths{
      distance_three_target,
      distance_four_target,
      fs::path(distance_three_manifest_path),
      fs::path(distance_four_manifest_path),
  };
  for (const auto& path : initial_output_paths) {
    if (!path.parent_path().empty()) {
      fs::create_directories(path.parent_path());
    }
  }

  const auto selected_manifest_shell = [&]() {
    Manifest selected;
    selected.format_version = 2;
    selected.length_limit = manifest.length_limit;
    selected.protocol_length_filter = manifest.protocol_length_filter;
    selected.base_support_only = false;
    selected.sources = manifest.sources;
    return selected;
  };
  auto distance_three_manifest = selected_manifest_shell();
  auto distance_four_manifest = selected_manifest_shell();
  std::map<std::uint32_t, std::uint64_t> distance_three_histogram;
  std::map<std::uint32_t, std::uint64_t> distance_four_histogram;

  const auto build_started = Clock::now();
  {
    std::ofstream distance_three_output(
        distance_three_temporary, std::ios::binary | std::ios::trunc);
    std::ofstream distance_four_output(
        distance_four_temporary, std::ios::binary | std::ios::trunc);
    if (!distance_three_output || !distance_four_output) {
      throw std::runtime_error("could not open temporary initial label caches");
    }
    const auto write_header = [](std::ofstream& output, std::uint32_t distance) {
      output.write(kLabelCacheMagic.data(), kLabelCacheMagic.size());
      write_little(output, std::uint32_t{1});
      write_little(output, distance);
      write_little(output, std::uint64_t{0});
      write_little(output, std::uint64_t{0});
    };
    write_header(distance_three_output, 3);
    write_header(distance_four_output, 4);

    struct InitialBases {
      std::vector<utsp::Mask> distance_three;
      std::vector<utsp::Mask> distance_four;
    };
    constexpr std::size_t minimum_chunk = 4096;
    const auto chunk_size = std::max<std::size_t>(
        minimum_chunk, static_cast<std::size_t>(requested_workers) * 1024U);
    for (std::size_t chunk_begin = 0;
         chunk_begin < manifest.supports.size();
         chunk_begin += chunk_size) {
      const auto chunk_end =
          std::min(manifest.supports.size(), chunk_begin + chunk_size);
      const auto count = chunk_end - chunk_begin;
      std::vector<InitialBases> bases(count);
      std::atomic<std::size_t> next{0};
      const auto worker_count = std::min<std::size_t>(requested_workers, count);
      std::vector<std::exception_ptr> errors(worker_count);
      std::vector<std::thread> threads;
      threads.reserve(worker_count);
      for (std::size_t worker = 0; worker < worker_count; ++worker) {
        threads.emplace_back([&, worker] {
          while (true) {
            const auto offset = next.fetch_add(1);
            if (offset >= count) return;
            try {
              const auto& record = manifest.supports[chunk_begin + offset];
              auto& output = bases[offset];
              if (source_schur) {
                if (record.source_index >= source_caches.size()) {
                  throw std::invalid_argument(
                      "pointed record source index is out of range");
                }
                output.distance_three = source_schur_pointed_quotient_basis(
                    source_caches[record.source_index], record);
              } else {
                output.distance_three =
                    utsp::build_logical_label_quotient_basis(
                        record.points,
                        record.ambient_dimension,
                        false,
                        3);
              }
              output.distance_four =
                  utsp::restrict_logical_label_quotient_basis_minimum_distance(
                      record.points,
                      record.ambient_dimension,
                      output.distance_three,
                      4,
                      false);
            } catch (...) {
              errors[worker] = std::current_exception();
              next.store(count);
              return;
            }
          }
        });
      }
      for (auto& thread : threads) thread.join();
      for (const auto& error : errors) {
        if (error) std::rethrow_exception(error);
      }
      for (std::size_t offset = 0; offset < count; ++offset) {
        const auto& record = manifest.supports[chunk_begin + offset];
        const auto& distance_three = bases[offset].distance_three;
        const auto& distance_four = bases[offset].distance_four;
        if (distance_three.size() > 64 || distance_four.size() > 64) {
          throw std::logic_error("initial label quotient dimension exceeds 64");
        }
        ++distance_three_histogram[
            static_cast<std::uint32_t>(distance_three.size())];
        ++distance_four_histogram[
            static_cast<std::uint32_t>(distance_four.size())];
        if (distance_three.size() >= 5) {
          distance_three_manifest.supports.push_back(record);
          write_little(
              distance_three_output,
              static_cast<std::uint8_t>(distance_three.size()));
          for (const auto row : distance_three) {
            write_little(distance_three_output, row);
          }
        }
        if (!distance_four.empty()) {
          distance_four_manifest.supports.push_back(record);
          write_little(
              distance_four_output,
              static_cast<std::uint8_t>(distance_four.size()));
          for (const auto row : distance_four) {
            write_little(distance_four_output, row);
          }
        }
      }
    }

    const auto finalize_header = [](std::ofstream& output,
                                    const Manifest& selected) {
      output.seekp(16, std::ios::beg);
      write_little(output, static_cast<std::uint64_t>(selected.supports.size()));
      write_little(output, manifest_cache_fingerprint(selected));
      output.seekp(0, std::ios::end);
      output.flush();
      if (!output) {
        throw std::runtime_error("failed to flush temporary initial label cache");
      }
    };
    finalize_header(distance_three_output, distance_three_manifest);
    finalize_header(distance_four_output, distance_four_manifest);
  }
  write_manifest(distance_three_manifest, distance_three_manifest_path);
  write_manifest(distance_four_manifest, distance_four_manifest_path);
  const auto promote = [](const fs::path& temporary, const fs::path& target) {
    std::error_code error;
    fs::rename(temporary, target, error);
    if (error) {
      fs::remove(temporary);
      throw std::runtime_error(
          "could not promote initial label cache: " + error.message());
    }
  };
  promote(distance_three_temporary, distance_three_target);
  promote(distance_four_temporary, distance_four_target);
  const auto build_finished = Clock::now();

  const auto print_profile = [](std::uint32_t distance,
                                std::uint32_t minimum_dimension,
                                const fs::path& cache,
                                const fs::path& selected_manifest,
                                const Manifest& selected,
                                const auto& histogram) {
    std::cout << "{\"minimum_distance\":" << distance
              << ",\"minimum_selected_quotient_dimension\":"
              << minimum_dimension
              << ",\"supports\":" << selected.supports.size()
              << ",\"cache_bytes\":" << fs::file_size(cache)
              << ",\"positive_manifest_bytes\":"
              << fs::file_size(selected_manifest)
              << ",\"quotient_dimension_histogram\":[";
    bool first = true;
    for (const auto& [dimension, count] : histogram) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << dimension << ',' << count << ']';
    }
    std::cout << "]}";
  };
  std::cout << "{\"schema\":\"utsp-native-initial-label-cache-summary-v1\"," 
            << "\"status\":\"complete\"," 
             << "\"input_supports\":" << manifest.supports.size() << ','
             << "\"workers\":" << requested_workers << ','
            << "\"source_schur\":"
            << (source_schur ? "true" : "false") << ','
            << "\"source_cache_seconds\":" << std::setprecision(12)
            << source_cache_seconds << ','
             << "\"distance_three\":";
  print_profile(
      3,
      5,
      distance_three_target,
      fs::path(distance_three_manifest_path),
      distance_three_manifest,
      distance_three_histogram);
  std::cout << ",\"distance_four\":";
  print_profile(
      4,
      1,
      distance_four_target,
      fs::path(distance_four_manifest_path),
      distance_four_manifest,
      distance_four_histogram);
  std::cout << ",\"read_seconds\":" << std::setprecision(12)
            << seconds(read_started, read_finished)
            << ",\"build_seconds\":"
            << seconds(build_started, build_finished) << "}\n";
  return 0;
}

int count_manifest(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto parent_text = argument_value(argc, argv, "--maximum-parent", false);
  const auto distance_text = argument_value(argc, argv, "--minimum-distance", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto work_chunk_text = argument_value(argc, argv, "--work-chunk", false);
  const auto status_seconds_text =
      argument_value(argc, argv, "--status-seconds", false);
  const auto start = start_text.empty() ? std::uint64_t{0} : parse_unsigned(start_text);
  const auto requested_count = count_text.empty()
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : parse_unsigned(count_text);
  const bool tensor_checksum = has_flag(argc, argv, "--tensor-checksum");
  const auto minimum_distance = distance_text.empty()
                                    ? std::uint32_t{3}
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(distance_text));
  const auto requested_workers = workers_text.empty()
                                     ? std::uint32_t{1}
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  const auto work_chunk = work_chunk_text.empty()
                              ? std::uint64_t{256}
                              : parse_unsigned(work_chunk_text);
  const bool trace_eligible = has_flag(argc, argv, "--trace-eligible");
  const auto status_seconds = status_seconds_text.empty()
                                  ? std::uint64_t{0}
                                  : parse_unsigned(status_seconds_text);
  if (requested_workers == 0) {
    throw std::invalid_argument("worker count must be positive");
  }
  if (work_chunk == 0) {
    throw std::invalid_argument("work chunk must be positive");
  }

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "manifest-count requires an origin-expanded or origin-reduced manifest");
  }
  const auto read_finished = Clock::now();
  const auto maximum_parent = parent_text.empty()
                                  ? manifest.length_limit
                                  : static_cast<std::uint32_t>(parse_unsigned(parent_text));
  if (start > manifest.supports.size()) {
    throw std::invalid_argument("manifest start is beyond the support count");
  }
  const auto available = static_cast<std::uint64_t>(manifest.supports.size()) - start;
  const auto selected_count = std::min(requested_count, available);

  struct CountTotals {
    std::map<std::uint32_t, std::uint64_t> label_dimensions;
    std::uint64_t isotropic_count = 0;
    std::uint64_t checksum = 0;
    std::uint64_t original_index_checksum = 0;
    double label_seconds = 0;
    double enumeration_seconds = 0;
    std::uint64_t eligible_supports = 0;
  };
  const auto worker_count = static_cast<std::uint32_t>(
      std::min<std::uint64_t>(requested_workers, selected_count));
  std::vector<CountTotals> worker_totals(worker_count);
  std::atomic<std::uint64_t> next_offset{0};
  std::mutex error_mutex;
  std::mutex trace_mutex;
  std::exception_ptr worker_error;
  const auto no_active_record = std::numeric_limits<std::uint64_t>::max();
  std::vector<std::atomic<std::uint64_t>> active_records(worker_count);
  std::vector<std::atomic<std::uint32_t>> active_sources(worker_count);
  std::vector<std::atomic<std::uint32_t>> active_dimensions(worker_count);
  for (std::uint32_t index = 0; index < worker_count; ++index) {
    active_records[index].store(no_active_record);
    active_sources[index].store(0);
    active_dimensions[index].store(0);
  }
  std::mutex status_mutex;
  std::condition_variable status_condition;
  bool status_finished = false;
  const bool parent_bound_applies =
      !manifest.protocol_length_filter || !parent_text.empty();
  auto worker = [&](std::uint32_t worker_index) {
    auto& totals = worker_totals[worker_index];
    try {
      while (true) {
        const auto chunk_start = next_offset.fetch_add(work_chunk);
        if (chunk_start >= selected_count) return;
        const auto chunk_end =
            std::min(selected_count, chunk_start + work_chunk);
        for (auto offset = chunk_start; offset < chunk_end; ++offset) {
          const auto record_index = start + offset;
          const auto& record = manifest.supports[record_index];
          totals.original_index_checksum ^=
              splitmix64(record.original_record_index);
          const auto label_started = Clock::now();
          const auto label_space = utsp::build_logical_label_space(
              record.points, record.ambient_dimension, false, minimum_distance);
          const auto label_finished = Clock::now();
          totals.label_seconds += seconds(label_started, label_finished);
          ++totals.label_dimensions[label_space.quotient_dimension()];
          if (label_space.quotient_dimension() < q ||
              (parent_bound_applies &&
               record.points.size() + q > maximum_parent)) {
            continue;
          }
          ++totals.eligible_supports;
          active_sources[worker_index].store(record.source_index);
          active_dimensions[worker_index].store(
              label_space.quotient_dimension());
          active_records[worker_index].store(record_index);
          if (trace_eligible) {
            const std::lock_guard lock(trace_mutex);
            std::cerr << "{\"event\":\"enumeration_start\",\"record\":"
                      << record_index << ",\"worker\":" << worker_index
                      << ",\"source\":" << record.source_index
                      << ",\"ambient_dimension\":"
                      << static_cast<unsigned>(record.ambient_dimension)
                      << ",\"label_dimension\":"
                      << label_space.quotient_dimension() << "}\n";
          }
          const auto enumeration_started = Clock::now();
          std::uint64_t record_isotropic_count = 0;
          if (tensor_checksum) {
            const auto summary =
                utsp::summarize_isotropic_subspaces(label_space, q);
            record_isotropic_count = summary.subspace_count;
            totals.isotropic_count += record_isotropic_count;
            totals.checksum ^= splitmix64(
                summary.tensor_checksum ^ splitmix64(record_index) ^
                splitmix64(summary.subspace_count));
          } else {
            record_isotropic_count =
                utsp::enumerate_totally_isotropic_subspaces(label_space, q);
            totals.isotropic_count += record_isotropic_count;
          }
          const auto record_enumeration_seconds =
              seconds(enumeration_started, Clock::now());
          totals.enumeration_seconds += record_enumeration_seconds;
          if (trace_eligible) {
            const std::lock_guard lock(trace_mutex);
            std::cerr << "{\"event\":\"enumeration_finish\",\"record\":"
                      << record_index << ",\"worker\":" << worker_index
                      << ",\"isotropic_subspaces\":"
                      << record_isotropic_count << ",\"seconds\":"
                      << std::setprecision(12) << record_enumeration_seconds
                      << "}\n";
          }
          active_records[worker_index].store(no_active_record);
        }
      }
    } catch (...) {
      const std::lock_guard lock(error_mutex);
      if (!worker_error) worker_error = std::current_exception();
      next_offset.store(selected_count);
    }
  };
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  std::thread status_reporter;
  if (status_seconds != 0) {
    status_reporter = std::thread([&] {
      const auto interval = std::chrono::seconds(status_seconds);
      std::unique_lock status_lock(status_mutex);
      while (!status_condition.wait_for(
          status_lock, interval, [&] { return status_finished; })) {
        const auto assigned = std::min(
            selected_count, next_offset.load());
        status_lock.unlock();
        const std::lock_guard trace_lock(trace_mutex);
        std::cerr << "{\"event\":\"manifest_count_status\","
                  << "\"records_assigned\":" << assigned
                  << ",\"active\":[";
        bool first_active = true;
        for (std::uint32_t index = 0; index < worker_count; ++index) {
          const auto record_index = active_records[index].load();
          if (record_index == no_active_record) continue;
          if (!first_active) std::cerr << ',';
          first_active = false;
          std::cerr << "{\"worker\":" << index << ",\"record\":"
                    << record_index << ",\"source\":"
                    << active_sources[index].load()
                    << ",\"label_dimension\":"
                    << active_dimensions[index].load() << '}';
        }
        std::cerr << "]}\n";
        status_lock.lock();
      }
    });
  }
  for (std::uint32_t index = 0; index < worker_count; ++index) {
    workers.emplace_back(worker, index);
  }
  for (auto& thread : workers) thread.join();
  {
    const std::lock_guard lock(status_mutex);
    status_finished = true;
  }
  status_condition.notify_all();
  if (status_reporter.joinable()) status_reporter.join();
  if (worker_error) std::rethrow_exception(worker_error);

  std::map<std::uint32_t, std::uint64_t> label_dimensions;
  std::uint64_t isotropic_count = 0;
  std::uint64_t checksum = 0;
  std::uint64_t original_index_checksum = 0;
  double label_seconds = 0;
  double enumeration_seconds = 0;
  std::uint64_t eligible_supports = 0;
  for (const auto& totals : worker_totals) {
    for (const auto& [dimension, count] : totals.label_dimensions) {
      label_dimensions[dimension] += count;
    }
    isotropic_count += totals.isotropic_count;
    checksum ^= totals.checksum;
    original_index_checksum ^= totals.original_index_checksum;
    label_seconds += totals.label_seconds;
    enumeration_seconds += totals.enumeration_seconds;
    eligible_supports += totals.eligible_supports;
  }

  std::cout << "{\"schema\":\"utsp-native-manifest-count-v1\","
            << "\"source_spaces\":" << manifest.sources.size() << ','
            << "\"manifest_supports\":" << manifest.supports.size() << ','
            << "\"record_start\":" << start << ','
            << "\"records_processed\":" << selected_count << ','
            << "\"workers\":" << worker_count << ','
            << "\"length_limit\":" << maximum_parent << ','
            << "\"manifest_length_filter\":\""
            << (manifest.protocol_length_filter ? "protocol" : "parent")
            << "\","
            << "\"logical_dimension\":" << q << ','
            << "\"minimum_distance\":" << minimum_distance << ','
            << "\"eligible_supports\":" << eligible_supports << ','
            << "\"label_spaces_by_dimension\":[";
  bool first = true;
  for (const auto& [dimension, count] : label_dimensions) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << dimension << ',' << count << ']';
  }
  std::cout << "],\"isotropic_subspace_count\":" << isotropic_count << ','
            << "\"original_index_checksum\":\"" << std::hex
            << original_index_checksum << std::dec << "\","
            << "\"tensor_checksum_enabled\":"
            << (tensor_checksum ? "true" : "false") << ','
            << "\"tensor_checksum\":\"" << std::hex << checksum << std::dec
            << "\",\"read_seconds\":" << std::setprecision(12)
            << seconds(read_started, read_finished) << ','
            << "\"label_space_seconds\":" << label_seconds << ','
            << "\"enumeration_seconds\":" << enumeration_seconds << "}\n";
  return 0;
}

int catalogue_manifest(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto label_cache_path =
      argument_value(argc, argv, "--label-space-cache", false);
  const auto producer_binary_sha256 =
      argument_value(argc, argv, "--producer-binary-sha256", false);
  if (!producer_binary_sha256.empty() &&
      (producer_binary_sha256.size() != 64 ||
       !std::all_of(
           producer_binary_sha256.begin(),
           producer_binary_sha256.end(),
           [](unsigned char character) { return std::isxdigit(character); }))) {
    throw std::invalid_argument(
        "producer binary SHA-256 must contain 64 hexadecimal digits");
  }
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto distance_text = argument_value(argc, argv, "--minimum-distance", false);
  const auto dp_text = argument_value(argc, argv, "--maximum-dp-dimension", false);
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto support_workers_text =
      argument_value(argc, argv, "--support-workers", false);
  const auto zero_isotropic_workers_text =
      argument_value(argc, argv, "--zero-isotropic-workers", false);
  const auto raw_maximum_dimension_text =
      argument_value(argc, argv, "--raw-max-quotient-dimension", false);
  const auto minimum_quotient_dimension_text =
      argument_value(argc, argv, "--minimum-quotient-dimension", false);
  const auto maximum_quotient_dimension_text =
      argument_value(argc, argv, "--maximum-quotient-dimension", false);
  const auto start = start_text.empty() ? std::uint64_t{0}
                                        : parse_unsigned(start_text);
  const auto requested_count =
      count_text.empty() ? std::numeric_limits<std::uint64_t>::max()
                         : parse_unsigned(count_text);
  const auto workers = workers_text.empty()
                           ? std::uint32_t{1}
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  const auto support_workers = support_workers_text.empty()
                                   ? std::uint32_t{1}
                                   : static_cast<std::uint32_t>(
                                         parse_unsigned(
                                             support_workers_text));
  const bool require_zero_isotropic =
      has_flag(argc, argv, "--require-zero-isotropic");
  const auto zero_isotropic_workers =
      require_zero_isotropic
          ? (zero_isotropic_workers_text.empty()
                 ? std::uint32_t{1}
                 : static_cast<std::uint32_t>(
                       parse_unsigned(zero_isotropic_workers_text)))
          : std::uint32_t{0};
  if (!require_zero_isotropic && !zero_isotropic_workers_text.empty()) {
    throw std::invalid_argument(
        "--zero-isotropic-workers requires --require-zero-isotropic");
  }
  if (require_zero_isotropic && zero_isotropic_workers == 0) {
    throw std::invalid_argument(
        "zero-isotropic worker count must be positive");
  }
  const bool marked_code_orbits = has_flag(argc, argv, "--marked-orbits");
  const bool final_nondegenerate_orbits =
      has_flag(argc, argv, "--final-nondegenerate-orbits");
  const bool q3_predecessor_target_orbits =
      has_flag(argc, argv, "--q3-predecessor-target-orbits");
  const bool q3_predecessor_existence_orbits =
      has_flag(argc, argv, "--q3-predecessor-existence-orbits");
  const bool q5_q3_chain_cover_orbits =
      has_flag(argc, argv, "--q5-q3-chain-cover-orbits");
  const bool q5_q4_hitting_set_orbits =
      has_flag(argc, argv, "--q5-q4-hitting-set-orbits");
  const bool q7_primitive_orbits =
      has_flag(argc, argv, "--q7-primitive-orbits");
  const bool hybrid_orbits = has_flag(argc, argv, "--hybrid-orbits");
  const bool emit_positive_supports =
      has_flag(argc, argv, "--emit-positive-supports");
  if (static_cast<unsigned int>(marked_code_orbits) +
          static_cast<unsigned int>(final_nondegenerate_orbits) +
          static_cast<unsigned int>(q3_predecessor_target_orbits) +
          static_cast<unsigned int>(q3_predecessor_existence_orbits) +
          static_cast<unsigned int>(q5_q3_chain_cover_orbits) +
          static_cast<unsigned int>(q5_q4_hitting_set_orbits) +
          static_cast<unsigned int>(q7_primitive_orbits) +
          static_cast<unsigned int>(hybrid_orbits) >
      1U) {
    throw std::invalid_argument(
        "choose one marked-code enumeration mode");
  }
  const auto enumeration_mode =
      marked_code_orbits
          ? utsp::ProtocolEnumerationMode::marked_code_orbits
      : final_nondegenerate_orbits
          ? utsp::ProtocolEnumerationMode::
                final_dimension_nondegenerate_marked_code_orbits
      : q3_predecessor_target_orbits
          ? utsp::ProtocolEnumerationMode::
                q3_predecessor_target_marked_code_orbits
      : q3_predecessor_existence_orbits
          ? utsp::ProtocolEnumerationMode::q3_predecessor_target_existence
      : q5_q3_chain_cover_orbits
          ? utsp::ProtocolEnumerationMode::
                q5_q3_chain_cover_marked_code_orbits
      : q5_q4_hitting_set_orbits
          ? utsp::ProtocolEnumerationMode::
                q5_q4_hitting_set_marked_code_orbits
      : q7_primitive_orbits
          ? utsp::ProtocolEnumerationMode::q7_primitive_marked_code_orbits
      : hybrid_orbits ? utsp::ProtocolEnumerationMode::hybrid
                      : utsp::ProtocolEnumerationMode::raw_subspaces;
  const auto raw_maximum_quotient_dimension =
      raw_maximum_dimension_text.empty()
          ? std::uint32_t{14}
          : static_cast<std::uint32_t>(
                parse_unsigned(raw_maximum_dimension_text));
  const auto minimum_quotient_dimension =
      minimum_quotient_dimension_text.empty()
          ? std::uint32_t{0}
          : static_cast<std::uint32_t>(
                parse_unsigned(minimum_quotient_dimension_text));
  const auto maximum_quotient_dimension =
      maximum_quotient_dimension_text.empty()
          ? std::numeric_limits<std::uint32_t>::max()
          : static_cast<std::uint32_t>(
                parse_unsigned(maximum_quotient_dimension_text));
  const std::string_view enumeration_mode_name =
      marked_code_orbits ? "complete_marked_code_orbits"
      : final_nondegenerate_orbits
          ? "final_dimension_nondegenerate_marked_code_orbits"
      : q3_predecessor_target_orbits
          ? "q3_predecessor_target_marked_code_orbits"
      : q3_predecessor_existence_orbits
          ? "q3_predecessor_target_existence"
      : q5_q3_chain_cover_orbits
          ? "q5_q3_chain_cover_marked_code_orbits"
      : q5_q4_hitting_set_orbits
          ? "q5_q4_hitting_set_marked_code_orbits"
      : q7_primitive_orbits
          ? "complete_q7_primitive_marked_code_orbits"
      : hybrid_orbits    ? "hybrid_raw_and_marked_code_orbits"
                         : "raw_isotropic_subspaces";
  const auto minimum_distance = distance_text.empty()
                                    ? std::uint32_t{3}
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(distance_text));
  const auto maximum_dp_dimension = dp_text.empty()
                                        ? std::uint32_t{24}
                                        : static_cast<std::uint32_t>(
                                              parse_unsigned(dp_text));
  const auto read_started = Clock::now();
  auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "manifest-catalogue requires an origin-expanded or origin-reduced manifest");
  }
  if (start > manifest.supports.size()) {
    throw std::invalid_argument(
        "manifest catalogue start is beyond the support count");
  }
  const auto available =
      static_cast<std::uint64_t>(manifest.supports.size()) - start;
  const auto selected_count = std::min(requested_count, available);
  QuotientCacheSelection cached_quotient_bases;
  if (!label_cache_path.empty()) {
    cached_quotient_bases = read_quotient_cache(
        label_cache_path,
        manifest,
        start,
        selected_count,
        minimum_distance);
  }
  const auto read_finished = Clock::now();
  std::vector<utsp::ProtocolSupport> supports;
  supports.reserve(static_cast<std::size_t>(selected_count));
  for (std::uint64_t offset = 0; offset < selected_count; ++offset) {
    const auto index = start + offset;
    const auto& record = manifest.supports[index];
    const auto& source = manifest.sources[record.source_index];
    supports.push_back(utsp::ProtocolSupport{
        record.original_record_index,
        record.source_index,
        source.id,
        source.family,
        record.base_length,
        record.ambient_dimension,
        record.parity_case,
        record.origin,
        record.points,
    });
  }
  manifest.supports.clear();
  manifest.supports.shrink_to_fit();
  const auto search_started = Clock::now();
  const auto result = utsp::classify_protocol_frontier(
      supports,
      q,
      minimum_distance,
      maximum_dp_dimension,
      enumeration_mode,
      raw_maximum_quotient_dimension,
      workers,
      support_workers,
      zero_isotropic_workers,
      minimum_quotient_dimension,
      maximum_quotient_dimension,
      emit_positive_supports,
      cached_quotient_bases.view());
  const auto search_finished = Clock::now();

  namespace fs = std::filesystem;
  const fs::path target(output_path);
  if (!target.parent_path().empty()) fs::create_directories(target.parent_path());
  const fs::path temporary = target.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) {
      throw std::runtime_error("could not open protocol output: " +
                               temporary.string());
    }
    output << "{\n  \"schema\": \""
           << (emit_positive_supports
                   ? "utsp-native-protocol-frontier-shard-v3"
               : enumeration_mode !=
                             utsp::ProtocolEnumerationMode::raw_subspaces ||
                         require_zero_isotropic
                   ? "utsp-native-protocol-frontier-shard-v2"
                   : "utsp-native-protocol-frontier-shard-v1")
           << "\",\n"
           << "  \"status\": \"complete\",\n"
           << "  \"producer_binary_sha256\": ";
    if (producer_binary_sha256.empty()) {
      output << "null";
    } else {
      write_json_string(output, producer_binary_sha256);
    }
    output << ",\n"
           << "  \"matrix_scope\": \"full_projective\",\n"
           << "  \"enumeration_mode\": \"" << enumeration_mode_name
           << "\",\n"
           << "  \"raw_maximum_quotient_dimension\": "
           << raw_maximum_quotient_dimension << ",\n"
           << "  \"quotient_dimension_interval\": ["
           << minimum_quotient_dimension << ',';
    if (maximum_quotient_dimension ==
        std::numeric_limits<std::uint32_t>::max()) {
      output << "null";
    } else {
      output << maximum_quotient_dimension;
    }
    output << "],\n"
           << "  \"orbit_workers\": " << workers << ",\n"
           << "  \"support_workers\": " << support_workers << ",\n"
           << "  \"zero_isotropic_proof\": "
           << (require_zero_isotropic ? "true" : "false") << ",\n"
           << "  \"zero_isotropic_workers\": "
           << zero_isotropic_workers << ",\n"
           << "  \"input_manifest\": ";
    write_json_string(output, input_path);
    output << ",\n  \"label_space_cache_used\": "
           << (label_cache_path.empty() ? "false" : "true")
           << ",\n  \"manifest_length_filter\": \""
           << (manifest.protocol_length_filter ? "protocol" : "parent")
           << "\",\n  \"length_limit\": " << manifest.length_limit
           << ",\n  \"support_range\": {\"start\": " << start
           << ", \"count\": " << selected_count
           << ", \"end_exclusive\": " << (start + selected_count)
           << "}"
           << ",\n  \"logical_qubits\": " << q
           << ",\n  \"minimum_distance\": " << minimum_distance
           << ",\n  \"output_key_definition\": "
           << "\"lexicographically least interleaved symmetric-cubic tensor "
              "bits under GL(q,2), packed little-endian in 64-bit words\",\n"
           << "  \"radical_policy\": \"radical tensors are omitted from the "
              "magic frontier because deleting Clifford-trivial radical rows "
              "strictly improves S and cannot lower d_Z; they are "
           << (result.radical_statistics_exact
                   ? "structurally counted"
                   : "excluded on marked-orbit supports before enumeration")
           << "\",\n"
           << "  \"statistics\": {\n"
           << "    \"source_spaces\": " << manifest.sources.size() << ",\n"
           << "    \"supports_processed\": " << result.supports_processed
           << ",\n    \"quotient_dimension_filtered_supports\": "
           << result.quotient_dimension_filtered_supports
           << ",\n    \"eligible_supports\": " << result.eligible_supports
           << ",\n    \"raw_enumerated_supports\": "
           << result.raw_enumerated_supports
           << ",\n    \"marked_orbit_enumerated_supports\": "
           << result.marked_orbit_enumerated_supports
           << ",\n    \"isotropic_subspace_statistics_exact\": "
           << (result.isotropic_subspace_statistics_exact ? "true" : "false")
           << ",\n    \"isotropic_subspaces\": ";
    if (result.isotropic_subspace_statistics_exact) {
      output << result.isotropic_subspaces;
    } else {
      output << "null";
    }
    output << ",\n    \"radical_statistics_exact\": "
           << (result.radical_statistics_exact ? "true" : "false")
           << ",\n    \"nondegenerate_subspaces\": "
           << result.nondegenerate_subspaces
           << ",\n    \"marked_code_orbits\": "
           << result.marked_code_orbits
           << ",\n    \"exact_distance_calls\": "
           << result.exact_distance_calls
           << ",\n    \"raw_tensor_forms\": " << result.raw_tensor_forms
           << ",\n    \"canonical_output_orbits\": "
           << result.canonical_output_orbits
           << ",\n    \"radical_dimension_counts\": {";
    bool first = true;
    for (const auto& [dimension, count] : result.radical_dimension_counts) {
      if (!first) output << ',';
      first = false;
      output << "\n      \"" << dimension << "\": " << count;
    }
    if (!result.radical_dimension_counts.empty()) output << '\n';
    output << "    },\n    \"read_seconds\": " << std::setprecision(12)
           << seconds(read_started, read_finished)
           << ",\n    \"search_seconds\": "
           << seconds(search_started, search_finished) << "\n  }";

    if (emit_positive_supports) {
      output << ",\n  \"positive_supports\": [";
      first = true;
      for (const auto& positive : result.positive_supports) {
        if (positive.support_offset >= supports.size()) {
          throw std::logic_error("positive support offset is out of range");
        }
        if (!first) output << ',';
        first = false;
        output << "\n    {\"support_offset\": " << positive.support_offset
               << ", \"manifest_support_index\": "
               << (start + positive.support_offset)
               << ", \"original_record_index\": "
               << positive.record_index << ", \"output_keys\": [";
        bool key_first = true;
        for (const auto& key : positive.output_keys) {
          if (!key_first) output << ',';
          key_first = false;
          output << "{\"logical_qubits\": " << key.logical_qubits
                 << ", \"canonical_key_words\": [";
          for (std::size_t word = 0; word < key.words.size(); ++word) {
            if (word != 0) output << ',';
            write_json_string(output, fixed_hex(key.words[word]));
          }
          output << "]}";
        }
        output << "]}";
      }
      if (!result.positive_supports.empty()) output << '\n' << "  ";
      output << ']';
    }
    output << ",\n  \"pareto_protocols\": [";

    constexpr std::array<std::string_view, 5> parity_names{
        "even",
        "remove_origin",
        "add_origin",
        "hyperplane_even",
        "hyperplane_add_origin",
    };
    first = true;
    for (const auto& witness : result.pareto_witnesses) {
      if (!first) output << ',';
      first = false;
      output << "\n    {\n      \"output\": {\n"
             << "        \"logical_qubits\": "
             << witness.output.logical_qubits << ",\n"
             << "        \"canonical_key_words\": [";
      for (std::size_t word = 0; word < witness.output.words.size(); ++word) {
        if (word != 0) output << ',';
        output << '\n' << "          ";
        write_json_string(output, fixed_hex(witness.output.words[word]));
      }
      if (!witness.output.words.empty()) output << '\n' << "        ";
      output << "],\n        \"known_class_id\": ";
      if (witness.known_class_id.empty()) {
        output << "null";
      } else {
        write_json_string(output, witness.known_class_id);
      }
      output << "\n      },\n      \"d_Z\": " << witness.distance
             << ",\n      \"protocol_length_n\": "
             << witness.protocol_length
             << ",\n      \"space_footprint_S\": "
             << witness.space_footprint
             << ",\n      \"error_coefficient\": "
             << witness.error_coefficient
             << ",\n      \"source\": {\n        \"index\": "
             << witness.source_index << ",\n        \"id\": ";
      write_json_string(output, witness.source_id);
      output << ",\n        \"family\": ";
      write_json_string(output, witness.source_family);
      output << ",\n        \"space_length\": "
             << witness.source_space_length;
      output << ",\n        \"manifest_record_index\": "
             << witness.record_index
             << ",\n        \"ambient_dimension\": "
             << static_cast<std::uint32_t>(witness.ambient_dimension)
             << ",\n        \"parity_case\": ";
      write_json_string(output, parity_names.at(witness.parity_case));
      output << ",\n        \"origin\": ";
      if (witness.origin < 0) {
        output << "null";
      } else {
        output << witness.origin;
      }
      output << "\n      },\n      \"generator_rows\": [";
      bool row_first = true;
      for (const auto row : witness.canonical_logical_rows) {
        if (!row_first) output << ',';
        row_first = false;
        output << '\n' << "        ";
        write_json_string(output, row_bits(row, witness.points.size()));
      }
      for (std::uint32_t coordinate = 0;
           coordinate < witness.ambient_dimension; ++coordinate) {
        utsp::Mask row = 0;
        for (std::size_t column = 0; column < witness.points.size(); ++column) {
          if (((witness.points[column] >> coordinate) & 1U) != 0) {
            row |= utsp::Mask{1} << column;
          }
        }
        if (!row_first) output << ',';
        row_first = false;
        output << '\n' << "        ";
        write_json_string(output, row_bits(row, witness.points.size()));
      }
      if (!row_first) output << '\n' << "      ";
      output << "]\n    }";
    }
    if (!result.pareto_witnesses.empty()) output << '\n' << "  ";
    output << "]\n}\n";
    output.flush();
    if (!output) throw std::runtime_error("failed to write protocol output");
  }
  std::error_code error;
  fs::rename(temporary, target, error);
  if (error) {
    fs::remove(temporary);
    throw std::runtime_error("could not promote protocol output: " +
                             error.message());
  }

  std::cout << "{\"schema\":\"utsp-native-protocol-frontier-summary-v2\"," 
            << "\"input\":\"" << input_path << "\",\"output\":\""
            << output_path << "\",\"producer_binary_sha256\":";
  if (producer_binary_sha256.empty()) {
    std::cout << "null";
  } else {
    std::cout << '"' << producer_binary_sha256 << '"';
  }
  std::cout << ",\"logical_qubits\":" << q
            << ",\"enumeration_mode\":\"" << enumeration_mode_name
            << "\",\"label_space_cache_used\":"
            << (label_cache_path.empty() ? "false" : "true")
            << ",\"orbit_workers\":" << workers
            << ",\"support_workers\":" << support_workers
            << ",\"zero_isotropic_proof\":"
            << (require_zero_isotropic ? "true" : "false")
            << ",\"zero_isotropic_workers\":" << zero_isotropic_workers
            << ",\"support_start\":" << start
            << ",\"support_count\":" << selected_count
            << ",\"supports_processed\":" << result.supports_processed
            << ",\"quotient_dimension_filtered_supports\":"
            << result.quotient_dimension_filtered_supports
            << ",\"raw_enumerated_supports\":"
            << result.raw_enumerated_supports
            << ",\"marked_orbit_enumerated_supports\":"
            << result.marked_orbit_enumerated_supports
            << ",\"isotropic_subspaces\":";
  if (result.isotropic_subspace_statistics_exact) {
    std::cout << result.isotropic_subspaces;
  } else {
    std::cout << "null";
  }
  std::cout << ",\"nondegenerate_subspaces\":"
            << result.nondegenerate_subspaces
            << ",\"marked_code_orbits\":" << result.marked_code_orbits
            << ",\"raw_tensor_forms\":" << result.raw_tensor_forms
            << ",\"canonical_output_orbits\":"
            << result.canonical_output_orbits
            << ",\"positive_supports\":"
            << result.positive_supports.size()
            << ",\"pareto_protocols\":" << result.pareto_witnesses.size()
            << ",\"read_seconds\":"
            << seconds(read_started, read_finished)
            << ",\"search_seconds\":"
            << seconds(search_started, search_finished) << "}\n";
  return 0;
}

int radical_q8_manifest_check(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto distance_text =
      argument_value(argc, argv, "--minimum-distance", false);
  const auto start = start_text.empty() ? std::uint64_t{0}
                                        : parse_unsigned(start_text);
  const auto requested_count =
      count_text.empty() ? std::numeric_limits<std::uint64_t>::max()
                         : parse_unsigned(count_text);
  const auto workers = workers_text.empty()
                           ? std::uint32_t{1}
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  const auto minimum_distance = distance_text.empty()
                                    ? std::uint32_t{3}
                                    : static_cast<std::uint32_t>(
                                          parse_unsigned(distance_text));
  if (workers == 0) {
    throw std::invalid_argument("radical q=8 workers must be positive");
  }

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "radical q=8 check requires a pointed support manifest");
  }
  if (start > manifest.supports.size()) {
    throw std::invalid_argument(
        "radical q=8 start is beyond the support count");
  }
  const auto selected_count = std::min<std::uint64_t>(
      requested_count, manifest.supports.size() - start);
  const auto read_finished = Clock::now();

  std::uint64_t eligible = 0;
  std::uint64_t seed_orbits = 0;
  std::uint64_t parent_orbits = 0;
  std::uint64_t extension_vector_orbits = 0;
  std::uint64_t rejected_radical_dimension = 0;
  std::uint64_t duplicate_children = 0;
  std::uint64_t canonical_search_nodes = 0;
  std::vector<std::uint64_t> supports_with_witness;
  std::uint64_t witness_support =
      std::numeric_limits<std::uint64_t>::max();
  std::uint32_t witness_intrinsic_dimension = 0;
  std::vector<utsp::Mask> witness_quotient_rows;
  std::vector<utsp::Mask> witness_label_rows;
  std::array<std::uint64_t, 9> seeds_by_intrinsic{};
  const auto search_started = Clock::now();
  for (std::uint64_t offset = 0; offset < selected_count; ++offset) {
    const auto index = start + offset;
    const auto& record = manifest.supports[index];
    const auto label_space = utsp::build_logical_label_space(
        record.points,
        record.ambient_dimension,
        false,
        minimum_distance);
    if (label_space.quotient_dimension() < 8) continue;
    ++eligible;
    const auto result =
        utsp::search_q8_isotropic_subspaces_by_tensor_radical(
            label_space, workers);
    seed_orbits += result.seed_orbits;
    parent_orbits += result.parent_orbits;
    extension_vector_orbits += result.extension_vector_orbits;
    rejected_radical_dimension += result.rejected_radical_dimension;
    duplicate_children += result.duplicate_children;
    canonical_search_nodes += result.canonical_search_nodes;
    for (std::size_t dimension = 0;
         dimension < result.seed_orbits_by_intrinsic_dimension.size();
         ++dimension) {
      seeds_by_intrinsic[dimension] +=
          result.seed_orbits_by_intrinsic_dimension[dimension];
    }
    if (result.has_isotropic_subspace) {
      supports_with_witness.push_back(index);
      if (witness_quotient_rows.empty()) {
        if (result.witness_quotient_rows.size() != 8 ||
            utsp::gf2_rank(result.witness_quotient_rows) != 8) {
          throw std::logic_error(
              "radical q=8 search returned an invalid witness basis");
        }
        witness_support = index;
        witness_intrinsic_dimension =
            result.witness_intrinsic_dimension;
        witness_quotient_rows = result.witness_quotient_rows;
        witness_label_rows.reserve(8);
        for (const auto row : witness_quotient_rows) {
          witness_label_rows.push_back(
              utsp::linear_combination(row, label_space.quotient_basis));
        }
        if (utsp::tensor_radical_basis(witness_label_rows).size() !=
            8 - witness_intrinsic_dimension) {
          throw std::logic_error(
              "radical q=8 witness has the wrong tensor radical");
        }
        for (const auto& form : label_space.bilinear_form_rows) {
          for (const auto left : witness_quotient_rows) {
            const auto image = utsp::linear_combination(left, form);
            for (const auto right : witness_quotient_rows) {
              if ((std::popcount(image & right) & 1U) != 0) {
                throw std::logic_error(
                    "radical q=8 witness is not common-isotropic");
              }
            }
          }
        }
      }
    }
  }
  const auto search_finished = Clock::now();

  std::cout
      << "{\"schema\":\"utsp-q8-radical-stratified-existence-v1\","
      << "\"status\":\"complete\",\"input\":\"" << input_path
      << "\",\"support_start\":" << start
      << ",\"support_count\":" << selected_count
      << ",\"eligible_supports\":" << eligible
      << ",\"supports_with_q8_isotropic_subspace\":[";
  for (std::size_t index = 0; index < supports_with_witness.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << supports_with_witness[index];
  }
  std::cout << "],\"seed_orbits_by_intrinsic_dimension\":[";
  for (std::size_t dimension = 0; dimension < seeds_by_intrinsic.size();
       ++dimension) {
    if (dimension != 0) std::cout << ',';
    std::cout << seeds_by_intrinsic[dimension];
  }
  std::cout << "],\"witness_support\":";
  if (witness_quotient_rows.empty()) {
    std::cout << "null,\"witness_intrinsic_dimension\":null,"
              << "\"witness_quotient_rows\":[],"
              << "\"witness_label_rows\":[]";
  } else {
    std::cout << witness_support
              << ",\"witness_intrinsic_dimension\":"
              << witness_intrinsic_dimension
              << ",\"witness_quotient_rows\":";
    print_mask_array(witness_quotient_rows);
    std::cout << ",\"witness_label_rows\":";
    print_mask_array(witness_label_rows);
  }
  std::cout << ",\"seed_orbits\":" << seed_orbits
            << ",\"parent_orbits\":" << parent_orbits
            << ",\"extension_vector_orbits\":"
            << extension_vector_orbits
            << ",\"rejected_radical_dimension\":"
            << rejected_radical_dimension
            << ",\"duplicate_children\":" << duplicate_children
            << ",\"canonical_search_nodes\":"
            << canonical_search_nodes
            << ",\"read_seconds\":"
            << seconds(read_started, read_finished)
            << ",\"search_seconds\":"
            << seconds(search_started, search_finished) << "}\n";
  return 0;
}

int regress_manifest(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto parent_text = argument_value(argc, argv, "--maximum-parent", false);
  const auto start = start_text.empty() ? std::uint64_t{0} : parse_unsigned(start_text);
  const auto requested_count = count_text.empty()
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : parse_unsigned(count_text);
  const bool full_matrix = has_flag(argc, argv, "--full-matrix");
  const bool emit_metric_coverage =
      has_flag(argc, argv, "--emit-metric-coverage");

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "manifest-regression requires an origin-expanded or origin-reduced manifest");
  }
  const auto classifier_started = Clock::now();
  const utsp::KnownOutputClassifier classifier;
  const auto setup_finished = Clock::now();
  const auto maximum_parent = parent_text.empty()
                                  ? manifest.length_limit
                                  : static_cast<std::uint32_t>(parse_unsigned(parent_text));
  if (!full_matrix && manifest.protocol_length_filter && parent_text.empty()) {
    throw std::invalid_argument(
        "known-output check on a protocol-filtered manifest needs --maximum-parent");
  }
  if (q < 1 || q > 4) throw std::invalid_argument("regression q must be 1,...,4");
  if (start > manifest.supports.size()) {
    throw std::invalid_argument("manifest start is beyond the support count");
  }
  const auto available = static_cast<std::uint64_t>(manifest.supports.size()) - start;
  const auto selected_count = std::min(requested_count, available);

  std::map<std::uint32_t, std::uint64_t> label_dimensions;
  std::map<std::string_view, std::uint64_t> class_counts;
  std::map<std::uint32_t, std::uint64_t> distance_counts;
  std::uint64_t isotropic_count = 0;
  std::uint64_t intrinsic_count = 0;
  std::uint64_t feasible_count = 0;
  std::uint64_t exact_distance_calls = 0;
  std::uint64_t distance_checksum = 0;
  std::uint64_t eligible_supports = 0;
  using FrontierWitness = std::pair<std::uint32_t, std::uint64_t>;
  using FrontierPoint = std::pair<std::uint32_t, std::uint32_t>;
  using FrontierPartition = std::pair<std::string, std::uint32_t>;
  std::map<
      FrontierPartition,
      std::map<FrontierPoint, std::set<FrontierWitness>>>
      frontier;
  using CoverageMetric = std::tuple<
      std::string,
      std::uint32_t,
      std::uint32_t,
      std::uint32_t,
      std::uint64_t,
      std::uint32_t>;
  std::set<CoverageMetric> metric_coverage;
  double label_seconds = 0;
  double search_seconds = 0;

  for (std::uint64_t offset = 0; offset < selected_count; ++offset) {
    const auto& record = manifest.supports[start + offset];
    const auto label_started = Clock::now();
    const auto label_space = utsp::build_logical_label_space(
        record.points, record.ambient_dimension, false);
    label_seconds += seconds(label_started, Clock::now());
    ++label_dimensions[label_space.quotient_dimension()];
    if (label_space.quotient_dimension() < q ||
        (!full_matrix && record.points.size() + q > maximum_parent)) {
      continue;
    }
    ++eligible_supports;
    const bool contains_zero =
        std::find(record.points.begin(), record.points.end(), 0U) !=
        record.points.end();
    const utsp::DistanceContext distance_context(
        record.points, record.ambient_dimension);
    const auto search_started = Clock::now();
    isotropic_count += utsp::enumerate_totally_isotropic_subspaces(
        label_space, q,
        [&](std::span<const utsp::Mask> quotient_rows,
            std::span<const utsp::Mask> label_rows) {
          (void)quotient_rows;
          const auto signature = utsp::cubic_tensor_word(label_rows);
          const auto* output = classifier.classify(q, signature);
          if (output == nullptr || output->intrinsic_logical_qubits != q) return;
          ++intrinsic_count;
          ++class_counts[output->class_id];
          if (!full_matrix) {
            const auto completion = classifier.minimum_feasible_completion(
                q,
                signature,
                static_cast<std::uint32_t>(record.points.size()),
                contains_zero,
                maximum_parent);
            if (!completion) return;
          }
          ++feasible_count;
          ++exact_distance_calls;
          const auto distance = distance_context.evaluate(label_rows);
          if (!distance.distance || *distance.distance < 3) {
            throw std::logic_error("reverse candidate has distance below three");
          }
          ++distance_counts[*distance.distance];
          if (emit_metric_coverage) {
            metric_coverage.emplace(
                std::string(output->class_id),
                *distance.distance,
                static_cast<std::uint32_t>(record.points.size()),
                record.ambient_dimension + q,
                distance.error_coefficient,
                record.source_index);
          }
          const FrontierPartition partition{
              std::string(output->class_id), *distance.distance};
          const FrontierPoint point{
              static_cast<std::uint32_t>(record.points.size()),
              record.ambient_dimension + q};
          auto& points = frontier[partition];
          bool dominated = false;
          for (const auto& [existing, witnesses] : points) {
            (void)witnesses;
            if (existing.first <= point.first &&
                existing.second <= point.second && existing != point) {
              dominated = true;
              break;
            }
          }
          if (!dominated) {
            for (auto iterator = points.begin(); iterator != points.end();) {
              const auto& existing = iterator->first;
              if (point.first <= existing.first &&
                  point.second <= existing.second && existing != point) {
                iterator = points.erase(iterator);
              } else {
                ++iterator;
              }
            }
            points[point].insert(FrontierWitness{
                record.source_index, distance.error_coefficient});
          }
          distance_checksum ^= splitmix64(
              splitmix64(*distance.distance) ^
              splitmix64(distance.error_coefficient) ^
              splitmix64(start + offset) ^ splitmix64(signature));
        });
    search_seconds += seconds(search_started, Clock::now());
  }

  std::cout << "{\"schema\":\"utsp-native-known-output-regression-v1\","
            << "\"source_spaces\":" << manifest.sources.size() << ','
            << "\"manifest_supports\":" << manifest.supports.size() << ','
            << "\"record_start\":" << start << ','
            << "\"records_processed\":" << selected_count << ','
            << "\"maximum_parent_length\":" << maximum_parent << ','
            << "\"logical_dimension\":" << q << ','
            << "\"matrix_scope\":\""
            << (full_matrix ? "full_projective" : "same_ambient_completable")
            << "\"," 
            << "\"eligible_supports\":" << eligible_supports << ','
            << "\"label_spaces_by_dimension\":[";
  bool first = true;
  for (const auto& [dimension, count] : label_dimensions) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << dimension << ',' << count << ']';
  }
  std::cout << "],\"isotropic_subspaces\":" << isotropic_count
            << ",\"nonzero_intrinsic_outputs\":" << intrinsic_count
            << ",\"completion_feasible_outputs\":" << feasible_count
            << ",\"exact_distance_calls\":" << exact_distance_calls
            << ",\"class_counts\":{";
  first = true;
  for (const auto& [class_id, count] : class_counts) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '"' << class_id << "\":" << count;
  }
  std::cout << "},\"distance_counts\":{";
  first = true;
  for (const auto& [distance, count] : distance_counts) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '"' << distance << "\":" << count;
  }
  std::cout << "},\"pareto_metrics\":[";
  first = true;
  for (const auto& [partition, points] : frontier) {
    for (const auto& [point, witnesses] : points) {
      if (!first) std::cout << ',';
      first = false;
      std::set<std::uint64_t> coefficients;
      std::set<std::uint32_t> source_indices;
      for (const auto& [source_index, coefficient] : witnesses) {
        source_indices.insert(source_index);
        coefficients.insert(coefficient);
      }
      std::cout << "{\"class_id\":\"" << partition.first
                << "\",\"distance\":" << partition.second
                << ",\"protocol_length\":" << point.first
                << ",\"space_footprint\":" << point.second
                << ",\"error_coefficients\":[";
      bool metric_first = true;
      for (const auto coefficient : coefficients) {
        if (!metric_first) std::cout << ',';
        metric_first = false;
        std::cout << coefficient;
      }
      std::cout << "],\"source_indices\":[";
      metric_first = true;
      for (const auto source_index : source_indices) {
        if (!metric_first) std::cout << ',';
        metric_first = false;
        std::cout << source_index;
      }
      std::cout << "]}";
    }
  }
  std::cout << "],\"metric_coverage_count\":" << metric_coverage.size();
  if (emit_metric_coverage) {
    std::cout << ",\"metric_coverage\":[";
    bool coverage_first = true;
    for (const auto& metric : metric_coverage) {
      if (!coverage_first) std::cout << ',';
      coverage_first = false;
      std::cout << "[\"" << std::get<0>(metric) << "\"," 
                << std::get<1>(metric) << ',' << std::get<2>(metric) << ','
                << std::get<3>(metric) << ',' << std::get<4>(metric) << ','
                << std::get<5>(metric) << ']';
    }
    std::cout << ']';
  }
  std::cout << ",\"distance_checksum\":\"" << std::hex << distance_checksum
            << std::dec << "\",\"read_seconds\":" << std::setprecision(12)
            << seconds(read_started, classifier_started)
            << ",\"classifier_setup_seconds\":"
            << seconds(classifier_started, setup_finished)
            << ",\"label_space_seconds\":" << label_seconds
            << ",\"search_seconds\":" << search_seconds << "}\n";
  return 0;
}

int manifest_origin_gate(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto node_limit_text = argument_value(argc, argv, "--node-limit", false);
  const auto progress_text = argument_value(argc, argv, "--progress-every", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto start = start_text.empty() ? std::uint64_t{0} : parse_unsigned(start_text);
  const auto requested_count = count_text.empty()
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : parse_unsigned(count_text);
  const auto node_limit = node_limit_text.empty()
                              ? std::uint64_t{0}
                              : parse_unsigned(node_limit_text);
  const auto progress_every = progress_text.empty()
                                  ? std::uint64_t{0}
                                  : parse_unsigned(progress_text);
  const auto requested_workers = workers_text.empty()
                                     ? std::uint32_t{1}
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  if (requested_workers == 0) {
    throw std::invalid_argument("worker count must be positive");
  }
  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  const auto read_finished = Clock::now();
  if (start > manifest.sources.size()) {
    throw std::invalid_argument("source start exceeds manifest source count");
  }
  const auto selected_count = std::min<std::uint64_t>(
      requested_count, manifest.sources.size() - start);
  std::vector<const PointedRecord*> base_records(manifest.sources.size(), nullptr);
  for (const auto& record : manifest.supports) {
    if (record.parity_case == 0 && record.origin >= 0 &&
        base_records[record.source_index] == nullptr) {
      base_records[record.source_index] = &record;
    }
  }

  struct DimensionTotals {
    std::uint64_t sources = 0;
    std::uint64_t raw_origins = 0;
    std::uint64_t origin_orbits = 0;
    std::uint64_t invariant_buckets = 0;
    std::uint64_t exact_tests = 0;
    std::uint64_t exact_nodes = 0;
    std::uint64_t coarse_nodes = 0;
    std::uint64_t unmarked_nodes = 0;
    std::uint64_t refined_dual_nodes = 0;
    std::uint64_t primal_nodes = 0;
    std::uint64_t second_order_refinements = 0;
    std::uint64_t refined_sources = 0;
    std::uint64_t marked_relation_sources = 0;
    std::uint64_t primal_search_sources = 0;
  };
  std::map<std::uint32_t, DimensionTotals> by_dimension;
  struct SourceOutcome {
    std::uint32_t dimension = 0;
    std::uint64_t raw_origins = 0;
    std::uint64_t origin_orbits = 0;
    std::uint64_t invariant_buckets = 0;
    std::uint64_t exact_tests = 0;
    std::uint64_t exact_nodes = 0;
    std::uint64_t coarse_nodes = 0;
    std::uint64_t unmarked_nodes = 0;
    std::uint64_t refined_dual_nodes = 0;
    std::uint64_t primal_nodes = 0;
    std::uint64_t second_order_refinements = 0;
    bool refined = false;
    bool marked_relations = false;
    bool primal_search = false;
    bool complete = false;
    double elapsed_seconds = 0;
  };
  std::vector<SourceOutcome> outcomes(selected_count);
  std::atomic<std::uint64_t> next_source{0};
  std::atomic<std::uint64_t> finished_sources{0};
  std::mutex progress_mutex;
  std::mutex error_mutex;
  std::exception_ptr worker_error;
  const auto search_started = Clock::now();
  auto worker = [&]() {
    try {
      while (true) {
        const auto offset = next_source.fetch_add(1);
        if (offset >= selected_count) return;
        const auto source_index = static_cast<std::size_t>(start + offset);
        const auto* record = base_records[source_index];
        if (record == nullptr) {
          throw std::runtime_error("manifest source has no translated even support");
        }
        std::vector<std::uint32_t> base_points;
        base_points.reserve(record->points.size());
        for (const auto point : record->points) {
          base_points.push_back(
              point ^ static_cast<std::uint32_t>(record->origin));
        }
        std::sort(base_points.begin(), base_points.end());
        const auto source_started = Clock::now();
        const auto orbits = utsp::affine_origin_orbits(
            base_points, record->base_affine_dimension, node_limit);
        outcomes[offset] = SourceOutcome{
            record->base_affine_dimension,
            std::uint64_t{1} << record->base_affine_dimension,
            orbits.orbits.size(),
            orbits.invariant_bucket_count,
            orbits.exact_equivalence_tests,
            orbits.exact_search_nodes,
            orbits.coarse_search_nodes,
            orbits.unmarked_refinement_search_nodes,
            orbits.refined_dual_search_nodes,
            orbits.primal_search_nodes,
            orbits.second_order_refinements,
            orbits.used_refined_colors,
            orbits.used_marked_relations,
            orbits.used_primal_search,
            orbits.complete,
            seconds(source_started, Clock::now()),
        };
        const auto finished = finished_sources.fetch_add(1) + 1;
        if (progress_every != 0 && finished % progress_every == 0) {
          const std::lock_guard lock(progress_mutex);
          std::cerr << "origin-gate completed " << finished << '/'
                    << selected_count << " selected sources; latest="
                    << manifest.sources[source_index].id << '\n'
                    << std::flush;
        }
      }
    } catch (...) {
      const std::lock_guard lock(error_mutex);
      if (!worker_error) worker_error = std::current_exception();
      next_source.store(selected_count);
    }
  };
  const auto worker_count = static_cast<std::uint32_t>(
      std::min<std::uint64_t>(requested_workers, selected_count));
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (std::uint32_t index = 0; index < worker_count; ++index) {
    workers.emplace_back(worker);
  }
  for (auto& thread : workers) thread.join();
  if (worker_error) std::rethrow_exception(worker_error);
  const auto search_finished = Clock::now();

  bool complete = true;
  std::uint64_t incomplete_source_index = std::numeric_limits<std::uint64_t>::max();
  std::string incomplete_source_id;
  std::uint64_t completed_sources = 0;
  double maximum_source_seconds = 0;
  std::string maximum_source_id;
  for (std::size_t offset = 0; offset < outcomes.size(); ++offset) {
    const auto& outcome = outcomes[offset];
    const auto source_index = static_cast<std::size_t>(start + offset);
    auto& totals = by_dimension[outcome.dimension];
    ++totals.sources;
    totals.raw_origins += outcome.raw_origins;
    totals.origin_orbits += outcome.origin_orbits;
    totals.invariant_buckets += outcome.invariant_buckets;
    totals.exact_tests += outcome.exact_tests;
    totals.exact_nodes += outcome.exact_nodes;
    totals.coarse_nodes += outcome.coarse_nodes;
    totals.unmarked_nodes += outcome.unmarked_nodes;
    totals.refined_dual_nodes += outcome.refined_dual_nodes;
    totals.primal_nodes += outcome.primal_nodes;
    totals.second_order_refinements += outcome.second_order_refinements;
    totals.refined_sources += static_cast<std::uint64_t>(outcome.refined);
    totals.marked_relation_sources +=
        static_cast<std::uint64_t>(outcome.marked_relations);
    totals.primal_search_sources +=
        static_cast<std::uint64_t>(outcome.primal_search);
    if (outcome.elapsed_seconds > maximum_source_seconds) {
      maximum_source_seconds = outcome.elapsed_seconds;
      maximum_source_id = manifest.sources[source_index].id;
    }
    if (outcome.complete) {
      ++completed_sources;
    } else if (complete) {
      complete = false;
      incomplete_source_index = source_index;
      incomplete_source_id = manifest.sources[source_index].id;
    }
  }

  std::uint64_t raw_origins = 0;
  std::uint64_t origin_orbits = 0;
  std::uint64_t exact_tests = 0;
  std::uint64_t exact_nodes = 0;
  std::uint64_t coarse_nodes = 0;
  std::uint64_t unmarked_nodes = 0;
  std::uint64_t refined_dual_nodes = 0;
  std::uint64_t primal_nodes = 0;
  std::uint64_t second_order_refinements = 0;
  std::uint64_t refined_sources = 0;
  std::uint64_t marked_relation_sources = 0;
  std::uint64_t primal_search_sources = 0;
  std::cout << "{\"schema\":\"utsp-native-manifest-origin-gate-v1\"," 
            << "\"input\":\"" << input_path << "\",\"source_start\":"
            << start << ",\"sources_requested\":" << selected_count
            << ",\"sources_processed\":" << completed_sources
            << ",\"workers\":" << worker_count
            << ",\"dimensions\":[";
  bool first = true;
  for (const auto& [dimension, totals] : by_dimension) {
    if (!first) std::cout << ',';
    first = false;
    raw_origins += totals.raw_origins;
    origin_orbits += totals.origin_orbits;
    exact_tests += totals.exact_tests;
    exact_nodes += totals.exact_nodes;
    coarse_nodes += totals.coarse_nodes;
    unmarked_nodes += totals.unmarked_nodes;
    refined_dual_nodes += totals.refined_dual_nodes;
    primal_nodes += totals.primal_nodes;
    second_order_refinements += totals.second_order_refinements;
    refined_sources += totals.refined_sources;
    marked_relation_sources += totals.marked_relation_sources;
    primal_search_sources += totals.primal_search_sources;
    std::cout << "{\"ambient_dimension\":" << dimension
              << ",\"sources\":" << totals.sources
              << ",\"raw_origins\":" << totals.raw_origins
              << ",\"origin_orbits\":" << totals.origin_orbits
              << ",\"reduction_factor\":" << std::setprecision(12)
              << (totals.origin_orbits == 0
                      ? 0.0
                      : static_cast<double>(totals.raw_origins) /
                            totals.origin_orbits)
              << ",\"invariant_buckets\":" << totals.invariant_buckets
              << ",\"exact_equivalence_tests\":" << totals.exact_tests
              << ",\"exact_search_nodes\":" << totals.exact_nodes
              << ",\"coarse_search_nodes\":" << totals.coarse_nodes
              << ",\"unmarked_refinement_search_nodes\":"
              << totals.unmarked_nodes
              << ",\"refined_dual_search_nodes\":"
              << totals.refined_dual_nodes
              << ",\"primal_search_nodes\":" << totals.primal_nodes
              << ",\"second_order_refinements\":"
              << totals.second_order_refinements
              << ",\"refined_sources\":" << totals.refined_sources
              << ",\"marked_relation_sources\":"
              << totals.marked_relation_sources
              << ",\"primal_search_sources\":"
              << totals.primal_search_sources << '}';
  }
  std::cout << "],\"raw_origins\":" << raw_origins
            << ",\"origin_orbits\":" << origin_orbits
            << ",\"reduction_factor\":" << std::setprecision(12)
            << (origin_orbits == 0
                    ? 0.0
                    : static_cast<double>(raw_origins) / origin_orbits)
            << ",\"exact_equivalence_tests\":" << exact_tests
            << ",\"exact_search_nodes\":" << exact_nodes
            << ",\"coarse_search_nodes\":" << coarse_nodes
            << ",\"unmarked_refinement_search_nodes\":"
            << unmarked_nodes
            << ",\"refined_dual_search_nodes\":" << refined_dual_nodes
            << ",\"primal_search_nodes\":" << primal_nodes
            << ",\"second_order_refinements\":"
            << second_order_refinements
            << ",\"refined_sources\":" << refined_sources
            << ",\"marked_relation_sources\":"
            << marked_relation_sources
            << ",\"primal_search_sources\":" << primal_search_sources
            << ",\"maximum_source_seconds\":" << maximum_source_seconds
            << ",\"maximum_source_id\":\"" << maximum_source_id
            << "\",\"incomplete_source_index\":";
  if (complete) {
    std::cout << "null,\"incomplete_source_id\":null";
  } else {
    std::cout << incomplete_source_index << ",\"incomplete_source_id\":\""
              << incomplete_source_id << '\"';
  }
  std::cout
            << ",\"read_seconds\":" << seconds(read_started, read_finished)
            << ",\"search_seconds\":" << seconds(search_started, search_finished)
            << ",\"complete\":" << (complete ? "true" : "false")
            << "}\n";
  return complete ? 0 : 3;
}

[[nodiscard]] bool support_is_in_manifest_scope(
    const Manifest& manifest, std::size_t protocol_length) {
  const auto measured_length = protocol_length +
      static_cast<std::size_t>(!manifest.protocol_length_filter);
  return measured_length <= manifest.length_limit;
}

[[nodiscard]] std::size_t measured_manifest_length(
    const Manifest& manifest, std::size_t protocol_length) {
  return protocol_length +
      static_cast<std::size_t>(!manifest.protocol_length_filter);
}

[[nodiscard]] std::uint32_t point_rank(
    std::span<const std::uint32_t> points) {
  std::vector<utsp::Mask> vectors(points.begin(), points.end());
  return utsp::gf2_rank(vectors);
}

int representative_ledger_manifest(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto family = argument_value(argc, argv, "--family");
  const auto source_prefix = argument_value(argc, argv, "--source-prefix");
  const auto dimension = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--dimension")));
  const auto source_length = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--source-length")));
  const auto length_limit = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--length-limit")));
  const auto ledger_records =
      parse_unsigned(argument_value(argc, argv, "--ledger-records"));
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto global_start_text =
      argument_value(argc, argv, "--global-start", false);
  const auto width_text =
      argument_value(argc, argv, "--source-index-width", false);
  const auto start = start_text.empty() ? std::uint64_t{0}
                                        : parse_unsigned(start_text);
  const auto global_start = global_start_text.empty()
                                ? std::uint64_t{0}
                                : parse_unsigned(global_start_text);
  const auto source_index_width = width_text.empty()
                                      ? std::uint32_t{9}
                                      : static_cast<std::uint32_t>(
                                            parse_unsigned(width_text));
  if (dimension < 6 || dimension > 31) {
    throw std::invalid_argument(
        "representative-ledger dimension must be in 6..31");
  }
  if (source_length == 0 || source_length >= (std::uint32_t{1} << 8U) ||
      source_length > length_limit) {
    throw std::invalid_argument("invalid representative-ledger length scope");
  }
  if (source_index_width == 0 || source_index_width > 20) {
    throw std::invalid_argument("source index width must be in 1..20");
  }
  if (start > ledger_records) {
    throw std::invalid_argument("ledger start exceeds its record count");
  }
  const auto selected_count = count_text.empty()
                                  ? ledger_records - start
                                  : parse_unsigned(count_text);
  if (selected_count > ledger_records - start ||
      selected_count > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("selected ledger interval is invalid");
  }

  const auto mask_word_count = std::uint64_t{1} << (dimension - 6U);
  const auto record_word_count = std::uint64_t{5} + mask_word_count;
  if (record_word_count >
      std::numeric_limits<std::uint64_t>::max() / sizeof(std::uint64_t)) {
    throw std::invalid_argument("representative-ledger record is too wide");
  }
  const auto record_bytes = record_word_count * sizeof(std::uint64_t);
  const auto expected_bytes = ledger_records * record_bytes;
  if (ledger_records != 0 && expected_bytes / ledger_records != record_bytes) {
    throw std::invalid_argument("representative-ledger byte count overflows");
  }
  const auto actual_bytes = std::filesystem::file_size(input_path);
  if (actual_bytes != expected_bytes) {
    throw std::runtime_error(
        "representative-ledger byte count does not match its declared records");
  }

  const auto read_started = Clock::now();
  std::ifstream input(input_path, std::ios::binary);
  if (!input) {
    throw std::runtime_error(
        "could not open representative ledger: " + input_path);
  }
  if (start > static_cast<std::uint64_t>(
                  std::numeric_limits<std::streamoff>::max()) /
                  record_bytes) {
    throw std::invalid_argument("representative-ledger seek offset is too large");
  }
  input.seekg(static_cast<std::streamoff>(start * record_bytes));
  if (!input) throw std::runtime_error("could not seek representative ledger");

  Manifest output;
  output.length_limit = length_limit;
  output.protocol_length_filter = true;
  output.base_support_only = true;
  output.sources.reserve(static_cast<std::size_t>(selected_count));
  output.supports.reserve(static_cast<std::size_t>(selected_count));
  std::uint64_t member_count_sum = 0;
  std::uint64_t signature_checksum = 1469598103934665603ULL;
  std::uint64_t support_checksum = 1469598103934665603ULL;
  for (std::uint64_t local_index = 0; local_index < selected_count;
       ++local_index) {
    for (std::uint32_t word = 0; word < 4; ++word) {
      const auto value = read_little<std::uint64_t>(input);
      signature_checksum ^= splitmix64(value + word);
      signature_checksum *= 1099511628211ULL;
    }
    std::vector<std::uint32_t> points;
    points.reserve(source_length);
    for (std::uint64_t word_index = 0; word_index < mask_word_count;
         ++word_index) {
      auto word = read_little<std::uint64_t>(input);
      support_checksum ^= splitmix64(word + word_index);
      support_checksum *= 1099511628211ULL;
      while (word != 0) {
        const auto bit = static_cast<std::uint32_t>(std::countr_zero(word));
        points.push_back(
            static_cast<std::uint32_t>(word_index * 64U + bit));
        word &= word - 1U;
      }
    }
    const auto member_count = read_little<std::uint64_t>(input);
    if (member_count == 0) {
      throw std::runtime_error(
          "representative-ledger exact class has zero members");
    }
    if (std::numeric_limits<std::uint64_t>::max() - member_count_sum <
        member_count) {
      throw std::runtime_error("representative-ledger member sum overflows");
    }
    member_count_sum += member_count;
    if (points.size() != source_length) {
      throw std::runtime_error(
          "representative-ledger support has the wrong weight");
    }
    if (point_rank(points) != dimension) {
      throw std::runtime_error(
          "representative-ledger support lacks full affine rank");
    }
    utsp::validate_stabilizer_moments(points, dimension, true);

    const auto global_index = global_start + start + local_index;
    std::ostringstream id;
    id << source_prefix << std::setw(source_index_width) << std::setfill('0')
       << global_index + 1U;
    output.sources.push_back(SourceRecord{id.str(), family});
    output.supports.push_back(PointedRecord{
        static_cast<std::uint32_t>(local_index),
        static_cast<std::uint16_t>(source_length),
        static_cast<std::uint8_t>(dimension),
        static_cast<std::uint8_t>(dimension),
        0,
        0,
        std::move(points),
    });
  }
  const auto read_finished = Clock::now();
  const auto write_started = Clock::now();
  write_manifest(output, output_path);
  const auto write_finished = Clock::now();
  std::cout
      << "{\"schema\":\"utsp-native-representative-ledger-manifest-v1\"," 
      << "\"input\":\"" << input_path << "\",\"output\":\""
      << output_path << "\",\"dimension\":" << dimension
      << ",\"source_length\":" << source_length
      << ",\"length_limit\":" << length_limit
      << ",\"record_bytes\":" << record_bytes
      << ",\"ledger_records\":" << ledger_records
      << ",\"record_start\":" << start
      << ",\"records_converted\":" << selected_count
      << ",\"global_record_interval\":[" << global_start + start << ','
      << global_start + start + selected_count << ']'
      << ",\"member_count_sum\":" << member_count_sum
      << ",\"signature_checksum\":" << signature_checksum
      << ",\"support_checksum\":" << support_checksum
      << ",\"read_seconds\":" << seconds(read_started, read_finished)
      << ",\"write_seconds\":" << seconds(write_started, write_finished)
      << "}\n";
  return 0;
}

void emit_reduced_support(
    const Manifest& manifest,
    std::size_t minimum_measured_length,
    std::size_t maximum_measured_length,
    std::uint32_t source_index,
    std::uint16_t base_length,
    std::uint8_t base_affine_dimension,
    std::uint8_t ambient_dimension,
    std::uint8_t parity_case,
    std::int32_t origin,
    std::vector<std::uint32_t> points,
    std::set<std::pair<std::uint8_t, std::vector<std::uint32_t>>>& seen,
    std::vector<PointedRecord>& output) {
  std::sort(points.begin(), points.end());
  const auto measured_length = measured_manifest_length(manifest, points.size());
  if (!support_is_in_manifest_scope(manifest, points.size()) ||
      measured_length < minimum_measured_length ||
      measured_length > maximum_measured_length ||
      point_rank(points) != ambient_dimension) {
    return;
  }
  const auto key = std::pair{ambient_dimension, points};
  if (!seen.insert(key).second) return;
  output.push_back(PointedRecord{
      source_index,
      base_length,
      base_affine_dimension,
      ambient_dimension,
      parity_case,
      origin,
      std::move(points),
  });
}

[[nodiscard]] bool pointed_record_less(
    const PointedRecord& left,
    const PointedRecord& right) {
  return std::tie(
             left.source_index,
             left.ambient_dimension,
             left.points,
             left.parity_case,
             left.origin,
             left.base_length,
             left.base_affine_dimension) <
         std::tie(
             right.source_index,
             right.ambient_dimension,
             right.points,
             right.parity_case,
             right.origin,
             right.base_length,
             right.base_affine_dimension);
}

[[nodiscard]] bool pointed_record_equal(
    const PointedRecord& left,
    const PointedRecord& right) {
  return !pointed_record_less(left, right) && !pointed_record_less(right, left);
}

int reduce_manifest_origins(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto node_limit_text = argument_value(argc, argv, "--node-limit", false);
  const bool raw_origin_mode = has_flag(argc, argv, "--raw-origins");
  const auto progress_text = argument_value(argc, argv, "--progress-every", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto minimum_length_text =
      argument_value(argc, argv, "--minimum-protocol-length", false);
  const auto maximum_length_text =
      argument_value(argc, argv, "--maximum-protocol-length", false);
  const auto node_limit = node_limit_text.empty()
                              ? std::uint64_t{0}
                              : parse_unsigned(node_limit_text);
  const auto progress_every = progress_text.empty()
                                  ? std::uint64_t{0}
                                  : parse_unsigned(progress_text);
  const auto requested_workers = workers_text.empty()
                                     ? std::uint32_t{1}
                                     : static_cast<std::uint32_t>(
                                           parse_unsigned(workers_text));
  const auto start = start_text.empty() ? std::uint64_t{0} : parse_unsigned(start_text);
  const auto requested_count = count_text.empty()
                                   ? std::numeric_limits<std::uint64_t>::max()
                                   : parse_unsigned(count_text);
  if (requested_workers == 0) {
    throw std::invalid_argument("worker count must be positive");
  }
  if (raw_origin_mode && !node_limit_text.empty()) {
    throw std::invalid_argument(
        "--raw-origins and --node-limit are mutually exclusive");
  }

  const auto read_started = Clock::now();
  const auto input = read_manifest(input_path);
  const auto read_finished = Clock::now();
  if ((!minimum_length_text.empty() || !maximum_length_text.empty()) &&
      !input.protocol_length_filter) {
    throw std::invalid_argument(
        "protocol-length bounds require a protocol-filtered manifest");
  }
  const auto minimum_measured_length = minimum_length_text.empty()
                                           ? std::uint64_t{0}
                                           : parse_unsigned(minimum_length_text);
  const auto maximum_measured_length = maximum_length_text.empty()
                                           ? std::uint64_t{input.length_limit}
                                           : parse_unsigned(maximum_length_text);
  if (minimum_measured_length > maximum_measured_length ||
      maximum_measured_length > input.length_limit) {
    throw std::invalid_argument("invalid protocol-length interval");
  }
  if (start > input.sources.size()) {
    throw std::invalid_argument("source start exceeds manifest source count");
  }
  const auto selected_count = std::min<std::uint64_t>(
      requested_count, input.sources.size() - start);
  std::vector<const PointedRecord*> base_records(input.sources.size(), nullptr);
  for (const auto& record : input.supports) {
    if (record.parity_case == 0 && record.origin >= 0 &&
        base_records[record.source_index] == nullptr) {
      base_records[record.source_index] = &record;
    }
  }

  struct SourceReduction {
    std::uint32_t dimension = 0;
    std::uint64_t raw_origins = 0;
    std::uint64_t origin_orbits = 0;
    std::uint64_t exact_nodes = 0;
    std::uint64_t coarse_nodes = 0;
    std::uint64_t unmarked_nodes = 0;
    std::uint64_t refined_dual_nodes = 0;
    std::uint64_t primal_nodes = 0;
    std::uint64_t second_order_refinements = 0;
    bool refined = false;
    bool marked_relations = false;
    bool primal_search = false;
    bool raw_origin_fallback = false;
    std::vector<PointedRecord> records;
  };
  std::vector<SourceReduction> reductions(selected_count);
  std::atomic<std::uint64_t> next_source{0};
  std::atomic<std::uint64_t> finished_sources{0};
  std::mutex progress_mutex;
  std::mutex error_mutex;
  std::exception_ptr worker_error;
  const auto reduce_started = Clock::now();
  auto worker = [&]() {
    try {
      while (true) {
        const auto offset = next_source.fetch_add(1);
        if (offset >= selected_count) return;
        const auto source_index = static_cast<std::size_t>(start + offset);
        const auto* base_record = base_records[source_index];
        if (base_record == nullptr) {
          throw std::runtime_error(
              "manifest source has no translated even support: " +
              input.sources[source_index].id);
        }
        const auto dimension = base_record->base_affine_dimension;
        std::vector<std::uint32_t> base_points;
        base_points.reserve(base_record->points.size());
        for (const auto point : base_record->points) {
          base_points.push_back(
              point ^ static_cast<std::uint32_t>(base_record->origin));
        }
        std::sort(base_points.begin(), base_points.end());
        SourceReduction result;
        result.dimension = dimension;
        result.raw_origins = std::uint64_t{1} << dimension;
        std::vector<std::uint32_t> origin_representatives;
        if (raw_origin_mode) {
          result.raw_origin_fallback = true;
          result.origin_orbits = result.raw_origins;
        } else {
          const auto orbit_result =
              utsp::affine_origin_orbits(base_points, dimension, node_limit);
          result.raw_origin_fallback = !orbit_result.complete;
          result.origin_orbits = result.raw_origin_fallback
                                     ? result.raw_origins
                                     : orbit_result.orbits.size();
          result.exact_nodes = orbit_result.exact_search_nodes;
          result.coarse_nodes = orbit_result.coarse_search_nodes;
          result.unmarked_nodes =
              orbit_result.unmarked_refinement_search_nodes;
          result.refined_dual_nodes = orbit_result.refined_dual_search_nodes;
          result.primal_nodes = orbit_result.primal_search_nodes;
          result.second_order_refinements =
              orbit_result.second_order_refinements;
          result.refined = orbit_result.used_refined_colors;
          result.marked_relations = orbit_result.used_marked_relations;
          result.primal_search = orbit_result.used_primal_search;
          if (orbit_result.complete) {
            origin_representatives.reserve(orbit_result.orbits.size());
            for (const auto& orbit : orbit_result.orbits) {
              origin_representatives.push_back(orbit.representative);
            }
          }
        }
        std::set<std::pair<std::uint8_t, std::vector<std::uint32_t>>> seen;
        const auto emit_origin = [&](std::uint32_t origin) {
          std::vector<std::uint32_t> translated;
          translated.reserve(base_points.size());
          for (const auto point : base_points) {
            translated.push_back(point ^ origin);
          }
          emit_reduced_support(
              input,
              minimum_measured_length,
              maximum_measured_length,
              static_cast<std::uint32_t>(source_index),
              base_record->base_length,
              dimension,
              dimension,
              0,
              static_cast<std::int32_t>(origin),
              translated,
              seen,
              result.records);

          std::sort(translated.begin(), translated.end());
          const bool origin_in_support =
              std::binary_search(base_points.begin(), base_points.end(), origin);
          if (origin_in_support) {
            if (translated.empty() || translated.front() != 0) {
              throw std::logic_error("translated support lost its marked origin");
            }
            translated.erase(translated.begin());
          } else {
            translated.insert(translated.begin(), 0);
          }
          emit_reduced_support(
              input,
              minimum_measured_length,
              maximum_measured_length,
              static_cast<std::uint32_t>(source_index),
              base_record->base_length,
              dimension,
              dimension,
              static_cast<std::uint8_t>(origin_in_support ? 1 : 2),
              static_cast<std::int32_t>(origin),
              std::move(translated),
              seen,
              result.records);
        };
        if (result.raw_origin_fallback) {
          // All origins are an exact cover, either by explicit strategy or
          // when automorphism search hits its certified node bound.
          for (std::uint32_t origin = 0; origin < result.raw_origins; ++origin) {
            emit_origin(origin);
          }
        } else {
          for (const auto origin : origin_representatives) {
            emit_origin(origin);
          }
        }

        if (dimension >= 31) {
          throw std::runtime_error("hyperplane embedding exceeds uint32_t points");
        }
        const auto extension_bit = std::uint32_t{1} << dimension;
        std::vector<std::uint32_t> embedded;
        embedded.reserve(base_points.size());
        for (const auto point : base_points) {
          embedded.push_back(point | extension_bit);
        }
        emit_reduced_support(
            input,
            minimum_measured_length,
            maximum_measured_length,
            static_cast<std::uint32_t>(source_index),
            base_record->base_length,
            dimension,
            static_cast<std::uint8_t>(dimension + 1),
            3,
            -1,
            embedded,
            seen,
            result.records);
        embedded.push_back(0);
        emit_reduced_support(
            input,
            minimum_measured_length,
            maximum_measured_length,
            static_cast<std::uint32_t>(source_index),
            base_record->base_length,
            dimension,
            static_cast<std::uint8_t>(dimension + 1),
            4,
            -1,
            std::move(embedded),
            seen,
            result.records);
        std::sort(
            result.records.begin(),
            result.records.end(),
            [](const PointedRecord& left, const PointedRecord& right) {
              return pointed_record_less(left, right);
            });
        reductions[offset] = std::move(result);

        const auto finished = finished_sources.fetch_add(1) + 1;
        if (progress_every != 0 && finished % progress_every == 0) {
          const std::lock_guard lock(progress_mutex);
          std::cerr << "origin-reduce completed " << finished << '/'
                    << selected_count << " selected sources; latest="
                    << input.sources[source_index].id << '\n'
                    << std::flush;
        }
      }
    } catch (...) {
      const std::lock_guard lock(error_mutex);
      if (!worker_error) worker_error = std::current_exception();
      next_source.store(selected_count);
    }
  };

  const auto worker_count = static_cast<std::uint32_t>(
      std::min<std::uint64_t>(requested_workers, selected_count));
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  for (std::uint32_t index = 0; index < worker_count; ++index) {
    workers.emplace_back(worker);
  }
  for (auto& thread : workers) thread.join();
  if (worker_error) std::rethrow_exception(worker_error);
  const auto reduce_finished = Clock::now();

  Manifest output;
  output.length_limit = input.length_limit;
  output.protocol_length_filter = input.protocol_length_filter;
  output.base_support_only = false;
  output.sources = input.sources;
  std::array<std::uint64_t, 5> parity_counts{};
  std::map<std::uint32_t, std::array<std::uint64_t, 10>> dimensions;
  std::uint64_t raw_origins = 0;
  std::uint64_t origin_orbits = 0;
  std::uint64_t exact_nodes = 0;
  std::uint64_t coarse_nodes = 0;
  std::uint64_t unmarked_nodes = 0;
  std::uint64_t refined_dual_nodes = 0;
  std::uint64_t primal_nodes = 0;
  std::uint64_t second_order_refinements = 0;
  std::uint64_t refined_sources = 0;
  std::uint64_t marked_relation_sources = 0;
  std::uint64_t primal_search_sources = 0;
  std::uint64_t raw_origin_fallback_sources = 0;
  for (auto& reduction : reductions) {
    auto& totals = dimensions[reduction.dimension];
    ++totals[0];
    totals[1] += reduction.raw_origins;
    totals[2] += reduction.origin_orbits;
    totals[3] += reduction.second_order_refinements;
    totals[4] += reduction.unmarked_nodes;
    totals[5] += static_cast<std::uint64_t>(reduction.marked_relations);
    totals[6] += reduction.refined_dual_nodes;
    totals[7] += reduction.primal_nodes;
    totals[8] += static_cast<std::uint64_t>(reduction.primal_search);
    totals[9] += static_cast<std::uint64_t>(reduction.raw_origin_fallback);
    raw_origins += reduction.raw_origins;
    origin_orbits += reduction.origin_orbits;
    exact_nodes += reduction.exact_nodes;
    coarse_nodes += reduction.coarse_nodes;
    unmarked_nodes += reduction.unmarked_nodes;
    refined_dual_nodes += reduction.refined_dual_nodes;
    primal_nodes += reduction.primal_nodes;
    second_order_refinements += reduction.second_order_refinements;
    refined_sources += static_cast<std::uint64_t>(reduction.refined);
    marked_relation_sources +=
        static_cast<std::uint64_t>(reduction.marked_relations);
    primal_search_sources +=
        static_cast<std::uint64_t>(reduction.primal_search);
    raw_origin_fallback_sources +=
        static_cast<std::uint64_t>(reduction.raw_origin_fallback);
    for (auto& record : reduction.records) {
      ++parity_counts.at(record.parity_case);
      output.supports.push_back(std::move(record));
    }
  }
  const auto write_started = Clock::now();
  write_manifest(output, output_path);
  const auto write_finished = Clock::now();

  std::cout << "{\"schema\":\"utsp-native-origin-reduced-manifest-v1\"," 
            << "\"input\":\"" << input_path << "\",\"output\":\""
            << output_path << "\",\"source_spaces\":" << selected_count
            << ",\"total_source_spaces\":" << input.sources.size()
            << ",\"source_start\":" << start
            << ",\"minimum_protocol_length\":"
            << minimum_measured_length
            << ",\"maximum_protocol_length\":"
            << maximum_measured_length
            << ",\"input_manifest_kind\":\""
            << (input.base_support_only ? "base_supports" : "pointed_supports")
            << '"'
            << ",\"origin_strategy\":\""
            << (raw_origin_mode
                    ? "raw_origin_exact_cover"
                    : "automorphism_orbit_reduction_with_exact_fallback")
            << '"'
            << ",\"input_supports\":" << input.supports.size()
            << ",\"reduced_supports\":" << output.supports.size()
            << ",\"support_reduction_factor\":" << std::setprecision(12);
  if (input.base_support_only) {
    std::cout << "null";
  } else {
    std::cout << (output.supports.empty()
                      ? 0.0
                      : static_cast<double>(input.supports.size()) /
                            output.supports.size());
  }
  std::cout
            << ",\"raw_origins\":" << raw_origins
            << ",\"origin_orbits\":" << origin_orbits
            << ",\"origin_reduction_factor\":"
            << (origin_orbits == 0
                    ? 0.0
                    : static_cast<double>(raw_origins) / origin_orbits)
            << ",\"refined_sources\":" << refined_sources
            << ",\"exact_search_nodes\":" << exact_nodes
            << ",\"coarse_search_nodes\":" << coarse_nodes
            << ",\"unmarked_refinement_search_nodes\":"
            << unmarked_nodes
            << ",\"refined_dual_search_nodes\":" << refined_dual_nodes
            << ",\"primal_search_nodes\":" << primal_nodes
            << ",\"primal_search_sources\":" << primal_search_sources
            << ",\"raw_origin_fallback_sources\":"
            << raw_origin_fallback_sources
            << ",\"second_order_refinements\":"
            << second_order_refinements
            << ",\"marked_relation_sources\":"
            << marked_relation_sources
            << ",\"parity_counts\":[";
  for (std::size_t index = 0; index < parity_counts.size(); ++index) {
    if (index != 0) std::cout << ',';
    std::cout << parity_counts[index];
  }
  std::cout << "],\"dimensions\":[";
  bool first = true;
  for (const auto& [dimension, totals] : dimensions) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"ambient_dimension\":" << dimension
              << ",\"sources\":" << totals[0]
              << ",\"raw_origins\":" << totals[1]
              << ",\"origin_orbits\":" << totals[2]
              << ",\"second_order_refinements\":" << totals[3]
              << ",\"unmarked_refinement_search_nodes\":" << totals[4]
              << ",\"marked_relation_sources\":" << totals[5]
              << ",\"refined_dual_search_nodes\":" << totals[6]
              << ",\"primal_search_nodes\":" << totals[7]
              << ",\"primal_search_sources\":" << totals[8]
              << ",\"raw_origin_fallback_sources\":" << totals[9] << '}';
  }
  std::cout << "],\"read_seconds\":"
            << seconds(read_started, read_finished)
            << ",\"reduction_seconds\":"
            << seconds(reduce_started, reduce_finished)
            << ",\"write_seconds\":"
            << seconds(write_started, write_finished) << "}\n";
  return 0;
}

int filter_manifest_quotient_dimension(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto minimum_text =
      argument_value(argc, argv, "--minimum-quotient-dimension", false);
  const auto maximum_text =
      argument_value(argc, argv, "--maximum-quotient-dimension", false);
  const auto distance_text =
      argument_value(argc, argv, "--minimum-distance", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto minimum_dimension =
      minimum_text.empty()
          ? std::uint32_t{0}
          : static_cast<std::uint32_t>(parse_unsigned(minimum_text));
  const auto maximum_dimension =
      maximum_text.empty()
          ? std::numeric_limits<std::uint32_t>::max()
          : static_cast<std::uint32_t>(parse_unsigned(maximum_text));
  const auto minimum_distance =
      distance_text.empty()
          ? std::uint32_t{3}
          : static_cast<std::uint32_t>(parse_unsigned(distance_text));
  const auto requested_workers =
      workers_text.empty()
          ? std::uint32_t{1}
          : static_cast<std::uint32_t>(parse_unsigned(workers_text));
  if (minimum_dimension > maximum_dimension) {
    throw std::invalid_argument("invalid quotient-dimension interval");
  }
  if (requested_workers == 0) {
    throw std::invalid_argument("worker count must be positive");
  }

  const auto read_started = Clock::now();
  auto input = read_manifest(input_path);
  const auto read_finished = Clock::now();
  if (input.base_support_only) {
    throw std::invalid_argument(
        "quotient filtering requires an origin-expanded or reduced manifest");
  }

  const auto filter_started = Clock::now();
  std::vector<std::uint8_t> selected(input.supports.size(), 0);
  const auto worker_count = static_cast<std::size_t>(
      std::min<std::uint64_t>(requested_workers, input.supports.size()));
  std::vector<std::map<std::uint32_t, std::uint64_t>> histograms(
      worker_count);
  std::vector<std::exception_ptr> errors(worker_count);
  std::atomic<std::uint64_t> next{0};
  constexpr std::uint64_t work_chunk = 1024;
  std::vector<std::thread> threads;
  threads.reserve(worker_count);
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    threads.emplace_back([&, worker] {
      try {
        while (true) {
          const auto begin = next.fetch_add(work_chunk);
          if (begin >= input.supports.size()) return;
          const auto end = std::min<std::uint64_t>(
              input.supports.size(), begin + work_chunk);
          for (auto index = begin; index < end; ++index) {
            const auto& record = input.supports[index];
            const auto label_space = utsp::build_logical_label_space(
                record.points,
                record.ambient_dimension,
                false,
                minimum_distance);
            const auto dimension = label_space.quotient_dimension();
            ++histograms[worker][dimension];
            selected[index] = static_cast<std::uint8_t>(
                dimension >= minimum_dimension &&
                dimension <= maximum_dimension);
          }
        }
      } catch (...) {
        errors[worker] = std::current_exception();
        next.store(input.supports.size());
      }
    });
  }
  for (auto& thread : threads) thread.join();
  for (const auto& error : errors) {
    if (error) std::rethrow_exception(error);
  }
  const auto filter_finished = Clock::now();

  Manifest output;
  output.format_version = 2;
  output.length_limit = input.length_limit;
  output.protocol_length_filter = input.protocol_length_filter;
  output.base_support_only = false;
  output.sources = input.sources;
  output.supports.reserve(static_cast<std::size_t>(
      std::count(selected.begin(), selected.end(), std::uint8_t{1})));
  std::uint64_t original_index_checksum = 0;
  for (std::size_t index = 0; index < input.supports.size(); ++index) {
    if (!selected[index]) continue;
    auto record = std::move(input.supports[index]);
    if (record.original_record_index ==
        std::numeric_limits<std::uint64_t>::max()) {
      record.original_record_index = index;
    }
    original_index_checksum ^= splitmix64(record.original_record_index);
    output.supports.push_back(std::move(record));
  }
  std::map<std::uint32_t, std::uint64_t> histogram;
  for (const auto& local : histograms) {
    for (const auto& [dimension, count] : local) {
      histogram[dimension] += count;
    }
  }
  const auto write_started = Clock::now();
  write_manifest(output, output_path);
  const auto write_finished = Clock::now();

  std::cout
      << "{\"schema\":\"utsp-quotient-filtered-manifest-v1\","
      << "\"input\":\"" << input_path << "\",\"output\":\""
      << output_path << "\",\"input_format_version\":"
      << input.format_version << ",\"output_format_version\":2"
      << ",\"input_supports\":" << input.supports.size()
      << ",\"selected_supports\":" << output.supports.size()
      << ",\"minimum_quotient_dimension\":" << minimum_dimension
      << ",\"maximum_quotient_dimension\":";
  if (maximum_dimension == std::numeric_limits<std::uint32_t>::max()) {
    std::cout << "null";
  } else {
    std::cout << maximum_dimension;
  }
  std::cout << ",\"minimum_distance\":" << minimum_distance
            << ",\"workers\":" << worker_count
            << ",\"original_index_checksum\":\"" << std::hex
            << original_index_checksum << std::dec
            << "\",\"dimension_histogram\":[";
  bool first = true;
  for (const auto& [dimension, count] : histogram) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << '[' << dimension << ',' << count << ']';
  }
  std::cout << "],\"read_seconds\":"
            << seconds(read_started, read_finished)
            << ",\"filter_seconds\":"
            << seconds(filter_started, filter_finished)
            << ",\"write_seconds\":"
            << seconds(write_started, write_finished) << "}\n";
  return 0;
}

int canonicalize_manifest_supports(int argc, char** argv) {
  const auto input_path = argument_value(argc, argv, "--input");
  const auto output_path = argument_value(argc, argv, "--output");
  const auto start_text = argument_value(argc, argv, "--start", false);
  const auto count_text = argument_value(argc, argv, "--count", false);
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto start = start_text.empty() ? std::uint64_t{0}
                                        : parse_unsigned(start_text);
  const auto requested_count =
      count_text.empty() ? std::numeric_limits<std::uint64_t>::max()
                         : parse_unsigned(count_text);
  const auto requested_workers =
      workers_text.empty()
          ? std::uint32_t{1}
          : static_cast<std::uint32_t>(parse_unsigned(workers_text));
  if (requested_workers == 0) {
    throw std::invalid_argument("worker count must be positive");
  }
  if (!utsp::marked_code_canonicalizer_available()) {
    throw std::runtime_error(
        "support canonicalization requires a Bliss-enabled build");
  }

  const auto read_started = Clock::now();
  const auto manifest = read_manifest(input_path);
  const auto read_finished = Clock::now();
  if (manifest.base_support_only) {
    throw std::invalid_argument(
        "support canonicalization requires an origin-expanded manifest");
  }
  if (start >= manifest.supports.size()) {
    throw std::invalid_argument("support start lies outside the manifest");
  }
  const auto count = std::min<std::uint64_t>(
      requested_count, manifest.supports.size() - start);
  if (count == 0) {
    throw std::invalid_argument("support range must be nonempty");
  }

  struct CanonicalSupportRecord {
    std::uint64_t support_index = 0;
    std::uint64_t original_record_index = 0;
    std::uint32_t source_index = 0;
    std::uint32_t support_length = 0;
    std::uint32_t ambient_dimension = 0;
    std::vector<utsp::Mask> key_words;
    std::uint64_t automorphism_group_order = 0;
    std::uint64_t canonical_search_nodes = 0;
  };
  std::vector<CanonicalSupportRecord> records(count);
  std::vector<std::exception_ptr> errors;
  const auto worker_count = static_cast<std::size_t>(
      std::min<std::uint64_t>(requested_workers, count));
  errors.resize(worker_count);
  std::atomic<std::uint64_t> next{0};
  std::vector<std::thread> threads;
  threads.reserve(worker_count);
  const auto canonical_started = Clock::now();
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    threads.emplace_back([&, worker] {
      try {
        while (true) {
          const auto offset = next.fetch_add(1);
          if (offset >= count) return;
          const auto support_index = start + offset;
          const auto& support = manifest.supports[support_index];
          const auto label_space = utsp::build_logical_label_space(
              support.points, support.ambient_dimension, false, 3);
          auto canonical = utsp::canonicalize_marked_code(label_space, {});
          records[offset] = CanonicalSupportRecord{
              support_index,
              support.original_record_index,
              support.source_index,
              static_cast<std::uint32_t>(support.points.size()),
              support.ambient_dimension,
              std::move(canonical.key_words),
              canonical.automorphism_group_order,
              canonical.canonical_search_nodes};
        }
      } catch (...) {
        errors[worker] = std::current_exception();
        next.store(count);
      }
    });
  }
  for (auto& thread : threads) thread.join();
  for (const auto& error : errors) {
    if (error) std::rethrow_exception(error);
  }
  const auto canonical_finished = Clock::now();

  std::map<std::tuple<
               std::uint32_t,
               std::uint32_t,
               std::vector<utsp::Mask>>,
           std::uint64_t>
      multiplicities;
  std::uint64_t canonical_search_nodes = 0;
  for (const auto& record : records) {
    ++multiplicities[{
        record.support_length, record.ambient_dimension, record.key_words}];
    canonical_search_nodes += record.canonical_search_nodes;
  }

  namespace fs = std::filesystem;
  const fs::path destination(output_path);
  if (!destination.parent_path().empty()) {
    fs::create_directories(destination.parent_path());
  }
  const fs::path temporary = destination.string() + ".tmp";
  const auto write_started = Clock::now();
  {
    std::ofstream output(temporary, std::ios::trunc);
    if (!output) {
      throw std::runtime_error(
          "could not open temporary support-key output");
    }
    output
        << "{\n  \"schema\": \"utsp-marked-support-canonical-keys-v1\",\n"
        << "  \"status\": \"complete\",\n"
        << "  \"input\": \"" << input_path << "\",\n"
        << "  \"input_format_version\": " << manifest.format_version << ",\n"
        << "  \"support_range\": {\"start\": " << start
        << ", \"count\": " << count << ", \"end_exclusive\": "
        << (start + count) << "},\n"
        << "  \"records\": [\n";
    for (std::size_t index = 0; index < records.size(); ++index) {
      const auto& record = records[index];
      if (index != 0) output << ",\n";
      output << "    {\"support_index\": " << record.support_index
             << ", \"original_record_index\": "
             << record.original_record_index << ", \"source_index\": "
             << record.source_index << ", \"support_length\": "
             << record.support_length << ", \"ambient_dimension\": "
             << record.ambient_dimension << ", \"key_words\": [";
      for (std::size_t word = 0; word < record.key_words.size(); ++word) {
        if (word != 0) output << ',';
        output << '\"' << utsp::hex_mask(record.key_words[word]) << '\"';
      }
      output << "], \"automorphism_group_order\": "
             << record.automorphism_group_order
             << ", \"canonical_search_nodes\": "
             << record.canonical_search_nodes << '}';
    }
    output << "\n  ]\n}\n";
    output.flush();
    if (!output) {
      throw std::runtime_error("failed to write support canonical keys");
    }
  }
  std::error_code rename_error;
  fs::rename(temporary, destination, rename_error);
  if (rename_error) {
    fs::remove(temporary);
    throw std::runtime_error(
        "could not promote support canonical keys: " +
        rename_error.message());
  }
  const auto write_finished = Clock::now();

  std::uint64_t repeated_supports = 0;
  for (const auto& [key, multiplicity] : multiplicities) {
    (void)key;
    if (multiplicity > 1) repeated_supports += multiplicity;
  }
  std::cout
      << "{\"schema\":\"utsp-marked-support-canonical-key-summary-v1\""
      << ",\"input_supports\":" << manifest.supports.size()
      << ",\"supports_processed\":" << count
      << ",\"distinct_canonical_keys\":" << multiplicities.size()
      << ",\"supports_in_repeated_keys\":" << repeated_supports
      << ",\"canonical_search_nodes\":" << canonical_search_nodes
      << ",\"workers\":" << worker_count
      << ",\"read_seconds\":" << seconds(read_started, read_finished)
      << ",\"canonical_seconds\":"
      << seconds(canonical_started, canonical_finished)
      << ",\"write_seconds\":" << seconds(write_started, write_finished)
      << "}\n";
  return 0;
}

int merge_manifests(int argc, char** argv) {
  const auto inputs_text = argument_value(argc, argv, "--inputs");
  const auto output_path = argument_value(argc, argv, "--output");
  std::vector<std::string> input_paths;
  std::size_t start = 0;
  while (start <= inputs_text.size()) {
    const auto end = inputs_text.find(',', start);
    const auto path = inputs_text.substr(start, end - start);
    if (path.empty()) throw std::invalid_argument("manifest input path is empty");
    input_paths.push_back(path);
    if (end == std::string::npos) break;
    start = end + 1;
  }
  if (input_paths.empty()) throw std::invalid_argument("no manifests to merge");

  const auto started = Clock::now();
  Manifest output = read_manifest(input_paths.front());
  if (output.base_support_only) {
    throw std::invalid_argument("manifest-merge requires expanded/reduced inputs");
  }
  const auto same_sources = [](const auto& left, const auto& right) {
    if (left.size() != right.size()) return false;
    for (std::size_t index = 0; index < left.size(); ++index) {
      if (left[index].id != right[index].id ||
          left[index].family != right[index].family) {
        return false;
      }
    }
    return true;
  };
  for (std::size_t index = 1; index < input_paths.size(); ++index) {
    auto input = read_manifest(input_paths[index]);
    if (input.base_support_only ||
        input.length_limit != output.length_limit ||
        input.protocol_length_filter != output.protocol_length_filter ||
        !same_sources(input.sources, output.sources)) {
      throw std::invalid_argument("manifest inputs have incompatible metadata");
    }
    output.supports.insert(
        output.supports.end(),
        std::make_move_iterator(input.supports.begin()),
        std::make_move_iterator(input.supports.end()));
  }
  std::sort(output.supports.begin(), output.supports.end(), pointed_record_less);
  if (std::adjacent_find(
          output.supports.begin(), output.supports.end(), pointed_record_equal) !=
      output.supports.end()) {
    throw std::invalid_argument("manifest inputs contain duplicate records");
  }
  write_manifest(output, output_path);
  std::cout << "{\"schema\":\"utsp-native-manifest-merge-v1\","
            << "\"status\":\"complete\",\"input_manifests\":"
            << input_paths.size() << ",\"source_spaces\":"
            << output.sources.size() << ",\"supports\":"
            << output.supports.size() << ",\"output\":\"" << output_path
            << "\",\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return 0;
}

std::string canonical_point_key(const utsp::CanonicalPointSet& canonical) {
  std::ostringstream key;
  key << canonical.ambient_dimension << ':' << canonical.linear_rank << ':';
  for (const auto point : canonical.points) key << point << ',';
  return key.str();
}

int sparse_completion_gate() {
  const utsp::KnownOutputClassifier classifier;
  constexpr std::array<std::uint32_t, 5> expected_orbits{0, 2, 4, 10, 23};
  bool okay = true;
  std::uint64_t tensors_checked = 0;
  std::uint64_t nonzero_checked = 0;
  std::cout << "{\"schema\":\"utsp-native-sparse-completion-gate-v1\","
            << "\"maximum_completion_weight\":7,\"dimensions\":[";
  for (std::uint32_t q = 1; q <= 4; ++q) {
    if (q != 1) std::cout << ',';
    utsp::SparseCompletionDecoder decoder(q);
    const std::uint64_t tensor_count =
        std::uint64_t{1} << decoder.tensor_bit_count();
    std::set<std::string> canonical_forms;
    std::map<std::string, std::uint16_t> q4_key_to_reference_orbit;
    std::map<std::uint32_t, std::uint64_t> weights;
    std::uint64_t canonical_mismatches = 0;
    std::uint64_t weight_mismatches = 0;
    for (std::uint64_t signature = 0; signature < tensor_count; ++signature) {
      const auto completion = decoder.decode(signature, 7);
      if (!completion) {
        okay = false;
        continue;
      }
      ++tensors_checked;
      ++weights[completion->weight()];
      canonical_forms.insert(canonical_point_key(completion->canonical_form));
      if (signature == 0) continue;
      ++nonzero_checked;
      const auto* output = classifier.classify(q, signature);
      if (output == nullptr || output->t_count != completion->weight()) {
        ++weight_mismatches;
        okay = false;
      }
      if (q == 4 && output != nullptr) {
        const auto key = canonical_point_key(completion->canonical_form);
        const auto [found, inserted] = q4_key_to_reference_orbit.emplace(
            key, output->four_qubit_canonical_mask);
        if (!inserted && found->second != output->four_qubit_canonical_mask) {
          ++canonical_mismatches;
          okay = false;
        }
      }
    }
    if (canonical_forms.size() != expected_orbits[q]) okay = false;
    std::cout << "{\"logical_dimension\":" << q
              << ",\"tensor_count\":" << tensor_count
              << ",\"canonical_orbit_count\":" << canonical_forms.size()
              << ",\"expected_canonical_orbit_count\":" << expected_orbits[q]
              << ",\"weight_mismatches\":" << weight_mismatches
              << ",\"q4_orbit_mapping_mismatches\":" << canonical_mismatches
              << ",\"weight_histogram\":{";
    bool first = true;
    for (const auto& [weight, count] : weights) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '"' << weight << "\":" << count;
    }
    std::cout << "}}";
  }
  std::cout << "],\"tensors_checked\":" << tensors_checked
            << ",\"nonzero_tensors_checked\":" << nonzero_checked
            << ",\"status\":\"" << (okay ? "pass" : "fail") << "\"}\n";
  return okay ? 0 : 1;
}

int sparse_decode(int argc, char** argv) {
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto signature = parse_unsigned(argument_value(argc, argv, "--signature"));
  const auto weight_text = argument_value(argc, argv, "--maximum-weight", false);
  const auto maximum_weight = weight_text.empty()
                                  ? std::uint32_t{7}
                                  : static_cast<std::uint32_t>(
                                        parse_unsigned(weight_text));
  utsp::SparseCompletionDecoder decoder(q);
  const auto completion = decoder.decode(signature, maximum_weight);
  std::cout << "{\"schema\":\"utsp-native-sparse-decode-v1\","
            << "\"logical_dimension\":" << q << ",\"tensor_signature\":\""
            << std::hex << signature << std::dec << "\",\"maximum_weight\":"
            << maximum_weight << ",\"found\":"
            << (completion ? "true" : "false");
  if (completion) {
    std::cout << ",\"weight\":" << completion->weight() << ",\"columns\":[";
    for (std::size_t index = 0; index < completion->columns.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << completion->columns[index];
    }
    std::cout << "],\"canonical_rank\":"
              << completion->canonical_form.linear_rank
              << ",\"canonical_points\":[";
    for (std::size_t index = 0;
         index < completion->canonical_form.points.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << completion->canonical_form.points[index];
    }
    std::cout << "],\"source_basis\":[";
    for (std::size_t index = 0;
         index < completion->canonical_form.source_basis.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << completion->canonical_form.source_basis[index];
    }
    std::cout << ']';
  }
  std::cout << "}\n";
  return 0;
}

int direct_tensor_canonicalize(int argc, char** argv) {
  const auto q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--q")));
  const auto signature = parse_unsigned(argument_value(argc, argv, "--signature"));
  const auto started = Clock::now();
  const auto canonical = utsp::canonicalize_cubic_tensor_direct(q, signature);
  const auto finished = Clock::now();
  std::cout << "{\"schema\":\"utsp-native-direct-tensor-canonical-v1\","
            << "\"logical_dimension\":" << q << ",\"tensor_signature\":\""
            << std::hex << signature << "\",\"canonical_key\":\""
            << canonical.canonical_key << "\",\"canonical_key_words\":[";
  for (std::size_t index = 0; index < canonical.canonical_key_words.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '"' << canonical.canonical_key_words[index] << '"';
  }
  std::cout << std::dec << "],\"basis\":[";
  for (std::size_t index = 0;
       index < canonical.new_basis_in_old_coordinates.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << canonical.new_basis_in_old_coordinates[index];
  }
  std::cout << "],\"complete_bases_evaluated\":"
            << canonical.complete_bases_evaluated
            << ",\"prefix_branches_pruned\":"
            << canonical.prefix_branches_pruned << ",\"elapsed_seconds\":"
            << std::setprecision(12) << seconds(started, finished) << "}\n";
  return 0;
}

int direct_tensor_canonicalize_labels(int argc, char** argv) {
  const auto rows = parse_masks(argument_value(argc, argv, "--labels"));
  const auto started = Clock::now();
  const auto canonical = utsp::canonicalize_cubic_tensor_from_labels(rows);
  const auto finished = Clock::now();
  const auto replay = utsp::interleaved_tensor_words_from_labels(
      rows, canonical.new_basis_in_old_coordinates);
  if (replay != canonical.canonical_key_words) {
    throw std::logic_error("label tensor transporter replay failed");
  }
  std::cout << "{\"schema\":\"utsp-native-direct-label-tensor-canonical-v1\","
            << "\"logical_dimension\":" << rows.size()
            << ",\"canonical_key_words\":[" << std::hex;
  for (std::size_t index = 0; index < canonical.canonical_key_words.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '"' << canonical.canonical_key_words[index] << '"';
  }
  std::cout << std::dec << "],\"basis\":[";
  for (std::size_t index = 0;
       index < canonical.new_basis_in_old_coordinates.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << canonical.new_basis_in_old_coordinates[index];
  }
  std::cout << "],\"complete_bases_evaluated\":"
            << canonical.complete_bases_evaluated
            << ",\"prefix_branches_pruned\":"
            << canonical.prefix_branches_pruned << ",\"elapsed_seconds\":"
            << std::setprecision(12) << seconds(started, finished) << "}\n";
  return 0;
}

int direct_tensor_gate() {
  const utsp::KnownOutputClassifier classifier;
  constexpr std::array<std::uint32_t, 5> expected_orbits{0, 2, 4, 10, 23};
  bool okay = true;
  std::uint64_t tensors_checked = 0;
  std::uint64_t transporter_mismatches = 0;
  std::uint64_t orbit_mapping_mismatches = 0;
  std::uint64_t complete_bases = 0;
  std::uint64_t pruned = 0;
  const auto started = Clock::now();
  std::cout << "{\"schema\":\"utsp-native-direct-tensor-gate-v1\","
            << "\"dimensions\":[";
  for (std::uint32_t q = 1; q <= 4; ++q) {
    if (q != 1) std::cout << ',';
    const auto bits = q + q * (q - 1) / 2 + q * (q - 1) * (q - 2) / 6;
    const auto count = std::uint64_t{1} << bits;
    std::set<std::uint64_t> keys;
    std::map<std::uint64_t, std::uint16_t> key_to_reference_orbit;
    std::uint64_t local_transporter_mismatches = 0;
    std::uint64_t local_orbit_mismatches = 0;
    for (std::uint64_t signature = 0; signature < count; ++signature) {
      const auto canonical =
          utsp::canonicalize_cubic_tensor_direct(q, signature);
      ++tensors_checked;
      complete_bases += canonical.complete_bases_evaluated;
      pruned += canonical.prefix_branches_pruned;
      keys.insert(canonical.canonical_key);
      if (utsp::interleaved_tensor_key(
              q, signature, canonical.new_basis_in_old_coordinates) !=
          canonical.canonical_key) {
        ++local_transporter_mismatches;
        ++transporter_mismatches;
        okay = false;
      }
      if (signature != 0) {
        const auto* output = classifier.classify(q, signature);
        const auto [found, inserted] = key_to_reference_orbit.emplace(
            canonical.canonical_key, output->four_qubit_canonical_mask);
        if (!inserted && found->second != output->four_qubit_canonical_mask) {
          ++local_orbit_mismatches;
          ++orbit_mapping_mismatches;
          okay = false;
        }
      }
    }
    if (keys.size() != expected_orbits[q]) okay = false;
    std::cout << "{\"logical_dimension\":" << q
              << ",\"tensor_count\":" << count
              << ",\"canonical_orbit_count\":" << keys.size()
              << ",\"expected_canonical_orbit_count\":" << expected_orbits[q]
              << ",\"transporter_mismatches\":"
              << local_transporter_mismatches
              << ",\"reference_orbit_mapping_mismatches\":"
              << local_orbit_mismatches << '}';
  }
  std::cout << "],\"tensors_checked\":" << tensors_checked
            << ",\"complete_bases_evaluated\":" << complete_bases
            << ",\"prefix_branches_pruned\":" << pruned
            << ",\"transporter_mismatches\":" << transporter_mismatches
            << ",\"reference_orbit_mapping_mismatches\":"
            << orbit_mapping_mismatches << ",\"elapsed_seconds\":"
            << std::setprecision(12) << seconds(started, Clock::now())
            << ",\"status\":\"" << (okay ? "pass" : "fail") << "\"}\n";
  return okay ? 0 : 1;
}

int five_qubit_tensor_orbit_census(int argc, char** argv) {
  const bool direct_check = !has_flag(argc, argv, "--skip-direct-check");
  const auto started = Clock::now();
  const auto census =
      utsp::enumerate_five_qubit_tensor_orbits(direct_check);
  const auto finished = Clock::now();
  const bool okay =
      census.tensor_count == (std::uint64_t{1} << 25) &&
      census.generated_linear_group_order == census.expected_linear_group_order &&
      (!direct_check ||
       (census.direct_canonical_checks == census.orbits.size() &&
        census.direct_canonical_mismatches == 0));

  std::cout
      << "{\"schema\":\"utsp-native-q5-full-tensor-orbit-census-v1\"," 
      << "\"logical_dimension\":5,\"tensor_bit_count\":25,"
      << "\"tensor_count\":" << census.tensor_count
      << ",\"equivalence_group\":\"GL(5,2)\","
      << "\"group_generators\":[[2,4,8,16,1],[3,2,4,8,16]],"
      << "\"generated_linear_group_order\":"
      << census.generated_linear_group_order
      << ",\"expected_linear_group_order\":"
      << census.expected_linear_group_order
      << ",\"generator_transitions\":" << census.generator_transitions
      << ",\"orbit_count\":" << census.orbits.size()
      << ",\"direct_canonicalizer_enabled\":"
      << (direct_check ? "true" : "false")
      << ",\"direct_canonical_checks\":"
      << census.direct_canonical_checks
      << ",\"direct_canonical_mismatches\":"
      << census.direct_canonical_mismatches
      << ",\"orbit_counts_by_radical_dimension\":{";
  for (std::size_t dimension = 0;
       dimension < census.orbit_counts_by_radical_dimension.size();
       ++dimension) {
    if (dimension) std::cout << ',';
    std::cout << '\"' << dimension << "\":"
              << census.orbit_counts_by_radical_dimension[dimension];
  }
  std::cout << "},\"tensor_counts_by_radical_dimension\":{";
  for (std::size_t dimension = 0;
       dimension < census.tensor_counts_by_radical_dimension.size();
       ++dimension) {
    if (dimension) std::cout << ',';
    std::cout << '\"' << dimension << "\":"
              << census.tensor_counts_by_radical_dimension[dimension];
  }
  std::cout << "},\"orbits\":[";
  for (std::size_t index = 0; index < census.orbits.size(); ++index) {
    if (index) std::cout << ',';
    const auto& orbit = census.orbits[index];
    std::cout << "{\"orbit_index\":" << index
              << ",\"canonical_key_words\":[\"0x" << std::hex
              << std::setw(16) << std::setfill('0') << orbit.canonical_key
              << "\"],\"canonical_standard_signature\":\"0x"
              << std::setw(16) << orbit.canonical_standard_signature
              << "\"" << std::setfill(' ') << std::dec
              << ",\"orbit_size\":" << orbit.orbit_size
              << ",\"stabilizer_size\":" << orbit.stabilizer_size
              << ",\"radical_dimension\":" << orbit.radical_dimension
              << '}';
  }
  std::cout << "],\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, finished) << ",\"status\":\""
            << (okay ? "pass" : "fail") << "\"}\n";
  return okay ? 0 : 1;
}

void write_tensor_extension_owner_map(
    const utsp::TensorExtensionCensus& census,
    const std::string& path) {
  namespace fs = std::filesystem;
  if (census.orbit_index_by_state.size() != census.extension_state_count) {
    throw std::logic_error("extension owner map was not retained");
  }
  const fs::path output_path(path);
  if (!output_path.parent_path().empty()) {
    fs::create_directories(output_path.parent_path());
  }
  const fs::path temporary = output_path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("could not open owner map output");
    constexpr std::array<char, 8> magic{
        'U', 'T', 'S', 'P', 'O', 'M', '1', '\0'};
    output.write(magic.data(), magic.size());
    write_little(output, std::uint32_t{1});
    write_little(output, census.parent_logical_qubits);
    write_little(output, census.parent_standard_signature);
    write_little(output, census.extension_bit_count);
    write_little(output, census.extension_state_count);
    write_little(
        output, static_cast<std::uint32_t>(census.orbits.size()));
    for (const auto owner : census.orbit_index_by_state) {
      write_little(output, owner);
    }
  }
  fs::rename(temporary, output_path);
}

void write_tensor_extension_transporter_map(
    const utsp::TensorExtensionCensus& census,
    const std::string& path) {
  namespace fs = std::filesystem;
  if (census.transport_to_representative_basis_by_state.size() !=
      census.extension_state_count) {
    throw std::logic_error("extension transporter map was not retained");
  }
  const fs::path output_path(path);
  if (!output_path.parent_path().empty()) {
    fs::create_directories(output_path.parent_path());
  }
  const fs::path temporary = output_path.string() + ".tmp";
  {
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
      throw std::runtime_error("could not open transporter map output");
    }
    constexpr std::array<char, 8> magic{
        'U', 'T', 'S', 'P', 'T', 'M', '1', '\0'};
    output.write(magic.data(), magic.size());
    write_little(output, std::uint32_t{1});
    write_little(output, census.parent_logical_qubits);
    write_little(output, census.parent_standard_signature);
    write_little(output, census.extension_bit_count);
    write_little(output, census.extension_state_count);
    write_little(
        output, static_cast<std::uint32_t>(census.orbits.size()));
    write_little(output, census.transporter_replay_checks);
    write_little(output, census.transporter_replay_mismatches);
    for (const auto basis :
         census.transport_to_representative_basis_by_state) {
      write_little(output, basis);
    }
  }
  fs::rename(temporary, output_path);
}

int tensor_extension_census(int argc, char** argv) {
  const auto parent_q = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--parent-q")));
  const auto parent_signature =
      parse_unsigned(argument_value(argc, argv, "--parent-signature"));
  const auto automorphism_order =
      parse_unsigned(argument_value(argc, argv, "--automorphism-order"));
  const auto generators = parse_basis_generators(
      argument_value(argc, argv, "--automorphisms", false), parent_q);
  const auto owner_map_path =
      argument_value(argc, argv, "--owner-map", false);
  const auto transporter_map_path =
      argument_value(argc, argv, "--transporter-map", false);
  const auto started = Clock::now();
  const auto census = utsp::census_tensor_extensions(
      parent_q,
      parent_signature,
      generators,
      automorphism_order,
      !owner_map_path.empty() || !transporter_map_path.empty(),
      !transporter_map_path.empty());
  if (!owner_map_path.empty()) {
    write_tensor_extension_owner_map(census, owner_map_path);
  }
  if (!transporter_map_path.empty()) {
    write_tensor_extension_transporter_map(census, transporter_map_path);
  }
  const auto finished = Clock::now();
  const auto orbit_mass = std::accumulate(
      census.orbits.begin(), census.orbits.end(), std::uint64_t{0},
      [](std::uint64_t total, const utsp::TensorExtensionOrbit& orbit) {
        return total + orbit.orbit_size;
      });
  const bool sizes_divide = std::all_of(
      census.orbits.begin(), census.orbits.end(),
      [&](const utsp::TensorExtensionOrbit& orbit) {
        return census.marked_stabilizer_group_order % orbit.orbit_size == 0;
      });
  const bool okay =
      orbit_mass == census.extension_state_count && sizes_divide &&
      census.actions.size() == generators.size() + parent_q &&
      census.transporter_replay_mismatches == 0 &&
      (owner_map_path.empty() ||
       census.orbit_index_by_state.size() == census.extension_state_count) &&
      (transporter_map_path.empty() ||
       census.transport_to_representative_basis_by_state.size() ==
           census.extension_state_count);

  std::cout
      << "{\"schema\":\"utsp-native-tensor-extension-census-v1\"," 
      << "\"status\":\"" << (okay ? "pass" : "fail") << "\"," 
      << "\"parent_logical_qubits\":" << parent_q
      << ",\"child_logical_qubits\":" << parent_q + 1
      << ",\"parent_standard_signature\":\"0x" << std::hex
      << parent_signature << "\"" << std::dec
      << ",\"extension_bit_count\":" << census.extension_bit_count
      << ",\"extension_state_count\":" << census.extension_state_count
      << ",\"parent_automorphism_group_order\":"
      << census.parent_automorphism_group_order
      << ",\"marked_stabilizer_group_order\":"
      << census.marked_stabilizer_group_order
      << ",\"parent_automorphism_generators\":[";
  for (std::size_t generator = 0; generator < generators.size(); ++generator) {
    if (generator) std::cout << ',';
    std::cout << '[';
    for (std::size_t coordinate = 0;
         coordinate < generators[generator].size(); ++coordinate) {
      if (coordinate) std::cout << ',';
      std::cout << generators[generator][coordinate];
    }
    std::cout << ']';
  }
  std::cout << "],\"induced_affine_actions\":[";
  for (std::size_t index = 0; index < census.actions.size(); ++index) {
    if (index) std::cout << ',';
    const auto& action = census.actions[index];
    std::cout << "{\"offset\":" << action.offset << ",\"columns\":[";
    for (std::size_t column = 0; column < action.columns.size(); ++column) {
      if (column) std::cout << ',';
      std::cout << action.columns[column];
    }
    std::cout << "]}";
  }
  std::cout << "],\"generator_transitions\":"
            << census.generator_transitions
            << ",\"owner_map_written\":"
            << (!owner_map_path.empty() ? "true" : "false")
            << ",\"transporter_map_written\":"
            << (!transporter_map_path.empty() ? "true" : "false")
            << ",\"transporter_replay_checks\":"
            << census.transporter_replay_checks
            << ",\"transporter_replay_mismatches\":"
            << census.transporter_replay_mismatches
            << ",\"marked_orbit_count\":" << census.orbits.size()
            << ",\"marked_orbit_mass\":" << orbit_mass
            << ",\"orbits\":[";
  for (std::size_t index = 0; index < census.orbits.size(); ++index) {
    if (index) std::cout << ',';
    const auto& orbit = census.orbits[index];
    const auto child_signature = utsp::combine_tensor_parent_and_extension(
        parent_q, parent_signature, orbit.representative);
    std::cout << "{\"extension_representative\":" << orbit.representative
              << ",\"child_standard_signature\":\"0x" << std::hex
              << child_signature << "\"" << std::dec
              << ",\"orbit_size\":" << orbit.orbit_size << '}';
  }
  std::cout << "],\"checks\":{" 
            << "\"every_extension_state_visited_once\":"
            << (orbit_mass == census.extension_state_count ? "true" : "false")
            << ",\"every_orbit_size_divides_marked_stabilizer_order\":"
            << (sizes_divide ? "true" : "false")
            << ",\"action_generator_count_is_parent_plus_translations\":"
            << (census.actions.size() == generators.size() + parent_q
                    ? "true"
                    : "false")
            << ",\"retained_owner_map_has_every_state\":"
            << (owner_map_path.empty() ||
                        census.orbit_index_by_state.size() ==
                            census.extension_state_count
                    ? "true"
                    : "false")
            << ",\"retained_transporters_replay_without_mismatch\":"
            << (transporter_map_path.empty() ||
                        (census.transporter_replay_checks != 0 &&
                         census.transporter_replay_mismatches == 0)
                    ? "true"
                    : "false")
            << "},\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, finished) << "}\n";
  return okay ? 0 : 1;
}

int tensor_unmark_q7(int argc, char** argv) {
  const auto authority = argument_value(argc, argv, "--authority");
  const auto marked_ledger = argument_value(argc, argv, "--marked-ledger");
  const auto q6_maps = argument_value(argc, argv, "--q6-map-directory");
  const auto q7_maps = argument_value(argc, argv, "--q7-owner-directory");
  const auto class_map = argument_value(argc, argv, "--class-map");
  const auto workers_text = argument_value(argc, argv, "--workers", false);
  const auto workers = workers_text.empty()
                           ? std::uint32_t{1}
                           : static_cast<std::uint32_t>(
                                 parse_unsigned(workers_text));
  const auto started = Clock::now();
  const auto result = utsp::unmark_q7_tensor_orbits(
      authority, marked_ledger, q6_maps, q7_maps, class_map, workers);
  std::cout
      << "{\"schema\":\"utsp-native-q7-tensor-unmarking-v1\"," 
      << "\"status\":\"" << (result.pass ? "pass" : "fail") << "\"," 
      << "\"q5_parent_count\":" << result.q5_parent_count
      << ",\"q6_parent_count\":" << result.q6_parent_count
      << ",\"q6_marked_orbit_count\":" << result.q6_marked_orbit_count
      << ",\"transitions_per_marked_orbit\":"
      << result.transitions_per_marked_orbit
      << ",\"q7_marked_orbit_count\":" << result.q7_marked_orbit_count
      << ",\"transition_edge_count\":" << result.transition_edge_count
      << ",\"expected_unmarked_orbit_count\":"
      << result.expected_unmarked_orbit_count
      << ",\"unmarked_orbit_count\":" << result.classes.size()
      << ",\"transition_affine_replay_checks\":"
      << result.transition_affine_replay_checks
      << ",\"transition_edge_replay_checks\":"
      << result.transition_edge_replay_checks
      << ",\"q7_tensor_mass\":" << result.q7_tensor_mass
      << ",\"gl7_order\":163849992929280"
      << ",\"class_map\":\"" << class_map << "\""
      << ",\"checks\":{"
      << "\"authority_and_marked_ledger_headers_match\":"
      << (result.checks.authority_and_marked_ledger_headers_match
              ? "true" : "false")
      << ",\"every_owner_and_transporter_map_header_matches\":"
      << (result.checks.every_owner_and_transporter_map_header_matches
              ? "true" : "false")
      << ",\"every_q6_marked_transport_replays\":"
      << (result.checks.every_q6_marked_transport_replays ? "true" : "false")
      << ",\"every_transition_affine_map_replays\":"
      << (result.checks.every_transition_affine_map_replays ? "true" : "false")
      << ",\"sampled_transition_edges_replay_as_full_tensors\":"
      << (result.checks.sampled_transition_edges_replay_as_full_tensors
              ? "true" : "false")
      << ",\"every_marked_parent_has_complete_extension_mass\":"
      << (result.checks.every_marked_parent_has_complete_extension_mass
              ? "true" : "false")
      << ",\"component_count_matches_independent_burnside_count\":"
      << (result.checks.component_count_matches_independent_burnside_count
              ? "true" : "false")
      << ",\"unmarked_orbit_mass_is_two_to_63\":"
      << (result.checks.unmarked_orbit_mass_is_two_to_63 ? "true" : "false")
      << ",\"every_orbit_size_divides_gl7_order\":"
      << (result.checks.every_orbit_size_divides_gl7_order ? "true" : "false")
      << ",\"class_map_written\":"
      << (result.checks.class_map_written ? "true" : "false")
      << "},\"classes\":[";
  for (std::size_t index = 0; index < result.classes.size(); ++index) {
    if (index) std::cout << ',';
    const auto& record = result.classes[index];
    std::cout << "{\"class_index\":" << record.class_index
              << ",\"representative_standard_signature\":\"0x"
              << std::hex << std::setw(16) << std::setfill('0')
              << record.representative_standard_signature << std::dec
              << std::setfill(' ') << "\",\"orbit_size\":"
              << record.orbit_size << ",\"stabilizer_size\":"
              << record.stabilizer_size << ",\"marked_node_count\":"
              << record.marked_node_count << '}';
  }
  std::cout << "],\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return result.pass ? 0 : 1;
}

int tensor_centroid(int argc, char** argv) {
  const auto rows = parse_masks(argument_value(argc, argv, "--labels"));
  const auto limit_text = argument_value(argc, argv, "--maximum-dimension", false);
  const auto limit = limit_text.empty()
                         ? std::uint32_t{20}
                         : static_cast<std::uint32_t>(parse_unsigned(limit_text));
  const auto started = Clock::now();
  const auto analysis = utsp::analyze_tensor_centroid(rows, limit);
  std::cout << "{\"schema\":\"utsp-native-tensor-centroid-v2\"," 
            << "\"logical_dimension\":" << rows.size()
            << ",\"centroid_dimension\":" << analysis.centroid_basis.size()
            << ",\"centroid_elements_tested\":"
            << analysis.centroid_elements_tested
            << ",\"enumeration_limit_exceeded\":"
            << (analysis.enumeration_limit_exceeded ? "true" : "false")
            << ",\"split_found\":" << (analysis.split ? "true" : "false");
  if (analysis.split) {
    std::cout << ",\"idempotent_words\":[" << std::hex;
    for (std::size_t index = 0;
         index < analysis.split->centroid_element.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << '\"' << analysis.split->centroid_element[index] << '\"';
    }
    std::cout << std::dec << "],\"image_dimension\":"
              << analysis.split->image_basis.size()
              << ",\"kernel_dimension\":"
              << analysis.split->kernel_basis.size()
              << ",\"image_basis\":[";
    for (std::size_t index = 0; index < analysis.split->image_basis.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << analysis.split->image_basis[index];
    }
    std::cout << "],\"kernel_basis\":[";
    for (std::size_t index = 0; index < analysis.split->kernel_basis.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << analysis.split->kernel_basis[index];
    }
    std::cout << ']';
  }
  std::cout << ",\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return 0;
}

int decomposed_tensor_canonicalize(int argc, char** argv) {
  const auto rows = parse_masks(argument_value(argc, argv, "--labels"));
  const auto started = Clock::now();
  const auto canonical = utsp::canonicalize_tensor_by_decomposition(rows);
  std::cout << "{\"schema\":\"utsp-native-decomposed-tensor-canonical-v1\","
            << "\"logical_dimension\":" << rows.size()
            << ",\"radical_dimension\":" << canonical.radical_basis.size()
            << ",\"component_count\":" << canonical.components.size()
            << ",\"components\":[";
  for (std::size_t component_index = 0;
       component_index < canonical.components.size(); ++component_index) {
    if (component_index) std::cout << ',';
    const auto& component = canonical.components[component_index];
    std::cout << "{\"dimension\":" << component.dimension
              << ",\"canonical_key_words\":[" << std::hex;
    for (std::size_t index = 0; index < component.canonical_key_words.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << '"' << component.canonical_key_words[index] << '"';
    }
    std::cout << std::dec << "],\"basis\":[";
    for (std::size_t index = 0;
         index < component.basis_in_original_coordinates.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << component.basis_in_original_coordinates[index];
    }
    std::cout << "]}";
  }
  std::cout << "],\"radical_basis\":[";
  for (std::size_t index = 0; index < canonical.radical_basis.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << canonical.radical_basis[index];
  }
  std::cout << "],\"canonical_basis\":[";
  for (std::size_t index = 0;
       index < canonical.canonical_basis_in_original_coordinates.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << canonical.canonical_basis_in_original_coordinates[index];
  }
  std::cout << "],\"centroid_elements_tested\":"
            << canonical.centroid_elements_tested << ",\"elapsed_seconds\":"
            << std::setprecision(12) << seconds(started, Clock::now()) << "}\n";
  return 0;
}

int hereditary_punctures(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto template_dimension = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--template-dimension")));
  const auto support = parse_points(argument_value(argc, argv, "--points"));
  const auto pattern = parse_points(argument_value(argc, argv, "--template"));
  const auto node_limit_text = argument_value(argc, argv, "--node-limit", false);
  const auto node_limit = node_limit_text.empty()
                              ? std::uint64_t{0}
                              : parse_unsigned(node_limit_text);
  const auto emit_text = argument_value(argc, argv, "--emit", false);
  const auto emit = emit_text.empty() ? std::uint64_t{0} : parse_unsigned(emit_text);
  const auto started = Clock::now();
  const auto result = utsp::find_hereditary_punctures(
      support, ambient, pattern, template_dimension, node_limit);
  std::cout << "{\"schema\":\"utsp-native-hereditary-punctures-v1\"," 
            << "\"ambient_dimension\":" << ambient
            << ",\"support_length\":" << support.size()
            << ",\"template_dimension\":" << template_dimension
            << ",\"template_size\":" << pattern.size()
            << ",\"complete\":" << (result.complete ? "true" : "false")
            << ",\"histogram_rejected\":"
            << (result.histogram_rejected ? "true" : "false")
            << ",\"search_nodes\":" << result.search_nodes
            << ",\"profile_rejections\":" << result.profile_rejections
            << ",\"complete_embeddings\":" << result.complete_embeddings
            << ",\"unique_direction_spaces\":"
            << result.unique_direction_spaces
            << ",\"unique_embedded_patterns\":"
            << result.unique_embedded_patterns
            << ",\"duplicate_embedded_patterns\":"
            << result.duplicate_embedded_patterns
            << ",\"candidate_cosets\":" << result.candidate_cosets
            << ",\"puncture_count\":" << result.puncture_index_masks.size()
            << ",\"punctures\":[";
  const auto emitted = std::min<std::uint64_t>(
      emit, result.puncture_index_masks.size());
  for (std::uint64_t puncture = 0; puncture < emitted; ++puncture) {
    if (puncture) std::cout << ',';
    const auto mask = result.puncture_index_masks[puncture];
    std::cout << "{\"mask\":\"" << utsp::hex_mask(mask) << "\",\"indices\":[";
    bool first = true;
    for (std::uint32_t index = 0; index < support.size(); ++index) {
      if (((mask >> index) & 1U) == 0) continue;
      if (!first) std::cout << ',';
      first = false;
      std::cout << index;
    }
    std::cout << "]}";
  }
  std::cout << "],\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return result.complete ? 0 : 3;
}

int origin_orbits(int argc, char** argv) {
  const auto ambient = static_cast<std::uint32_t>(
      parse_unsigned(argument_value(argc, argv, "--ambient")));
  const auto support = parse_points(argument_value(argc, argv, "--points"));
  const auto node_limit_text = argument_value(argc, argv, "--node-limit", false);
  const auto node_limit = node_limit_text.empty()
                              ? std::uint64_t{0}
                              : parse_unsigned(node_limit_text);
  const auto emit_text = argument_value(argc, argv, "--emit", false);
  const auto emit = emit_text.empty() ? std::uint64_t{0} : parse_unsigned(emit_text);
  const bool emit_transporters = has_flag(argc, argv, "--transporters");
  const auto started = Clock::now();
  const auto result = utsp::affine_origin_orbits(support, ambient, node_limit);
  std::cout << "{\"schema\":\"utsp-native-origin-orbits-v1\"," 
            << "\"ambient_dimension\":" << ambient
            << ",\"support_length\":" << support.size()
            << ",\"raw_origin_count\":" << (std::uint64_t{1} << ambient)
            << ",\"orbit_count\":" << result.orbits.size()
            << ",\"invariant_bucket_count\":"
            << result.invariant_bucket_count
            << ",\"exact_equivalence_tests\":"
            << result.exact_equivalence_tests
            << ",\"exact_search_nodes\":" << result.exact_search_nodes
            << ",\"coarse_search_nodes\":" << result.coarse_search_nodes
            << ",\"unmarked_refinement_search_nodes\":"
            << result.unmarked_refinement_search_nodes
            << ",\"refined_dual_search_nodes\":"
            << result.refined_dual_search_nodes
            << ",\"primal_search_nodes\":"
            << result.primal_search_nodes
            << ",\"used_refined_colors\":"
            << (result.used_refined_colors ? "true" : "false")
            << ",\"used_marked_relations\":"
            << (result.used_marked_relations ? "true" : "false")
            << ",\"used_primal_search\":"
            << (result.used_primal_search ? "true" : "false")
            << ",\"second_order_refinements\":"
            << result.second_order_refinements
            << ",\"complete\":" << (result.complete ? "true" : "false")
            << ",\"orbits\":[";
  const auto emitted = std::min<std::uint64_t>(emit, result.orbits.size());
  for (std::uint64_t index = 0; index < emitted; ++index) {
    if (index) std::cout << ',';
    const auto& orbit = result.orbits[index];
    std::cout << "{\"representative\":" << orbit.representative
              << ",\"size\":" << orbit.members.size()
              << ",\"members\":[";
    for (std::size_t member = 0; member < orbit.members.size(); ++member) {
      if (member) std::cout << ',';
      std::cout << orbit.members[member].origin;
    }
    std::cout << ']';
    if (emit_transporters) {
      std::cout << ",\"transporters\":[";
      for (std::size_t member = 0; member < orbit.members.size(); ++member) {
        if (member) std::cout << ',';
        const auto& transporter = orbit.members[member];
        std::cout << "{\"origin\":" << transporter.origin
                  << ",\"source_dual_basis\":[";
        for (std::size_t basis = 0;
             basis < transporter.source_dual_basis.size(); ++basis) {
          if (basis) std::cout << ',';
          std::cout << transporter.source_dual_basis[basis];
        }
        std::cout << "],\"image_dual_basis\":[";
        for (std::size_t basis = 0;
             basis < transporter.image_dual_basis.size(); ++basis) {
          if (basis) std::cout << ',';
          std::cout << transporter.image_dual_basis[basis];
        }
        std::cout << "]}";
      }
      std::cout << ']';
    }
    std::cout << '}';
  }
  std::cout << "],\"elapsed_seconds\":" << std::setprecision(12)
            << seconds(started, Clock::now()) << "}\n";
  return result.complete ? 0 : 3;
}

void usage() {
  std::cerr
      << "Usage:\n"
      << "  utsp-native self-test\n"
      << "  utsp-native inspect --ambient H --points p0,p1,... --q Q"
         " [--minimum-distance D] [--emit N]\n"
      << "  utsp-native support-automorphisms --ambient H"
         " --points p0,p1,...\n"
      << "  utsp-native support-vector-orbits --ambient H"
         " --points p0,p1,...\n"
      << "  utsp-native isotropic-orbit-cover --ambient H"
         " --points p0,p1,... --q Q [--orbit-depth D]\n"
      << "  utsp-native canonical-isotropic-orbits --ambient H"
         " --points p0,p1,... --q Q"
         " [--sector all|rank-at-most-one|even]"
         " [--target-signatures MASK,...] [--target-dimension Q]"
         " [--target-nondegenerate-seed K]"
         " [--complete-nondegenerate|--successor-nondegenerate"
         "|--primitive-nondegenerate"
         "|--q5-zero-hyperplane-primitive"
         "|--q5-zero-hyperplane-complete"
         "|--q5-q3-chain-cover"
         "|--q5-q3-chain-profile-presence"
         "|--q5-q4-hitting-set"
         "|--hereditary-nondegenerate"
         "|--retain-nondegenerate-from K]"
         " [--workers N] [--emit N] [--emit-orbit-data]\n"
      << "  utsp-native tensor-primitive-coverage --ambient H"
         " --points p0,p1,... --q Q"
         " [--hereditary-from K] [--workers N]\n"
      << "  utsp-native primitive-tensor-census --q Q"
         " [--sector all|zero|rank-one|full-alternating]"
         " [--emit N] [--direct-keys]\n"
      << "  utsp-native q7-primitive-tensor-authority [--workers N]\n"
      << "  utsp-native q8-primitive-tensor-authority-screen\n"
      << "  utsp-native manifest-count --input FILE --q Q [--maximum-parent C]"
         " [--minimum-distance D] [--start I] [--count N]"
         " [--workers N] [--work-chunk N] [--tensor-checksum]"
         " [--trace-eligible] [--status-seconds N]\n"
      << "  utsp-native manifest-label-cache --input FILE --output FILE"
         " [--positive-manifest FILE]"
         " [--minimum-selected-quotient-dimension W]"
         " [--minimum-distance D] [--workers N] [--source-schur]\n"
      << "  utsp-native manifest-initial-label-caches --input FILE"
         " --distance-three-cache FILE --distance-three-manifest FILE"
         " --distance-four-cache FILE --distance-four-manifest FILE"
         " [--workers N] [--source-schur]\n"
      << "  utsp-native manifest-label-cache-differential --input FILE"
         " --reference-cache FILE --candidate-cache FILE"
         " --minimum-distance D [--workers N]\n"
      << "  utsp-native manifest-catalogue --input FILE --output FILE --q Q"
         " [--producer-binary-sha256 HEX]"
         " [--label-space-cache FILE]"
         " [--minimum-distance D] [--maximum-dp-dimension S]"
         " [--start I] [--count N]"
         " [--marked-orbits|--final-nondegenerate-orbits|"
         "--q3-predecessor-target-orbits|"
         "--q3-predecessor-existence-orbits|"
         "--q5-q3-chain-cover-orbits|"
         "--q5-q4-hitting-set-orbits|"
         "--q7-primitive-orbits|--hybrid-orbits]"
         " [--raw-max-quotient-dimension W] [--workers N]"
          " [--support-workers N] [--minimum-quotient-dimension W]"
          " [--maximum-quotient-dimension W]"
          " [--emit-positive-supports]"
          " [--require-zero-isotropic --zero-isotropic-workers N]\n"
      << "  utsp-native manifest-q8-radical-check --input FILE"
         " [--start I] [--count N] [--workers N]"
         " [--minimum-distance D]\n"
      << "  utsp-native manifest-regression --input FILE --q Q"
         " [--maximum-parent C] [--start I] [--count N] [--full-matrix]"
         " [--emit-metric-coverage]\n"
      << "  utsp-native manifest-origin-gate --input FILE"
         " [--start I] [--count N] [--node-limit N]"
         " [--workers N] [--progress-every N]\n"
      << "  utsp-native manifest-origin-reduce --input FILE --output FILE"
         " [--start I] [--count N] [--node-limit N|--raw-origins]"
         " [--workers N]"
         " [--progress-every N] [--minimum-protocol-length N]"
         " [--maximum-protocol-length N]\n"
      << "  utsp-native representative-ledger-manifest --input FILE"
         " --output FILE --dimension M --source-length C --length-limit L"
         " --ledger-records N --family NAME --source-prefix PREFIX"
         " [--start I] [--count N] [--global-start I]"
         " [--source-index-width W]\n"
      << "  utsp-native manifest-filter-quotient --input FILE --output FILE"
         " [--minimum-quotient-dimension W]"
         " [--maximum-quotient-dimension W]"
         " [--minimum-distance D] [--workers N]\n"
      << "  utsp-native manifest-support-canonical-keys --input FILE"
         " --output FILE [--start I] [--count N] [--workers N]\n"
      << "  utsp-native manifest-merge --inputs FILE,... --output FILE\n"
      << "  utsp-native sparse-completion-gate\n"
      << "  utsp-native sparse-decode --q Q --signature MASK"
         " [--maximum-weight W]\n"
      << "  utsp-native direct-canonical --q Q --signature MASK\n"
      << "  utsp-native direct-canonical-labels --labels MASK,...\n"
      << "  utsp-native direct-tensor-gate\n"
      << "  utsp-native five-qubit-tensor-orbits [--skip-direct-check]\n"
      << "  utsp-native tensor-extension-census --parent-q Q"
         " --parent-signature MASK --automorphism-order N"
         " [--automorphisms b0,...;b0,...]"
         " [--owner-map FILE] [--transporter-map FILE]\n"
      << "  utsp-native tensor-unmark-q7 --authority FILE"
         " --marked-ledger FILE --q6-map-directory DIR"
         " --q7-owner-directory DIR --class-map FILE [--workers N]\n"
      << "  utsp-native tensor-centroid --labels MASK,..."
         " [--maximum-dimension C]\n"
      << "  utsp-native decomposed-canonical --labels MASK,...\n"
      << "  utsp-native hereditary-punctures --ambient M --points p0,..."
         " --template-dimension A --template r0,..."
         " [--node-limit N] [--emit N]\n"
      << "  utsp-native origin-orbits --ambient M --points p0,..."
         " [--node-limit N] [--emit N] [--transporters]\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      usage();
      return 2;
    }
    const std::string_view command = argv[1];
    if (command == "self-test") return self_test();
    if (command == "inspect") return inspect_support(argc, argv);
    if (command == "support-automorphisms") {
      return inspect_support_automorphisms(argc, argv);
    }
    if (command == "support-vector-orbits") {
      return inspect_support_vector_orbits(argc, argv);
    }
    if (command == "isotropic-orbit-cover") {
      return inspect_isotropic_orbit_cover(argc, argv);
    }
    if (command == "canonical-isotropic-orbits") {
      return inspect_canonical_isotropic_orbits(argc, argv);
    }
    if (command == "tensor-primitive-coverage") {
      return tensor_primitive_coverage(argc, argv);
    }
    if (command == "primitive-tensor-census") {
      return primitive_tensor_census(argc, argv);
    }
    if (command == "q7-primitive-tensor-authority") {
      return q7_primitive_tensor_authority(argc, argv);
    }
    if (command == "q8-primitive-tensor-authority-screen") {
      return q8_primitive_tensor_authority_screen();
    }
    if (command == "manifest-count") return count_manifest(argc, argv);
    if (command == "manifest-label-cache") {
      return cache_manifest_label_spaces(argc, argv);
    }
    if (command == "manifest-initial-label-caches") {
      return cache_manifest_initial_label_spaces(argc, argv);
    }
    if (command == "manifest-label-cache-differential") {
      return differential_label_space_caches(argc, argv);
    }
    if (command == "manifest-catalogue") {
      return catalogue_manifest(argc, argv);
    }
    if (command == "manifest-q8-radical-check") {
      return radical_q8_manifest_check(argc, argv);
    }
    if (command == "manifest-regression") return regress_manifest(argc, argv);
    if (command == "manifest-origin-gate") {
      return manifest_origin_gate(argc, argv);
    }
    if (command == "manifest-origin-reduce") {
      return reduce_manifest_origins(argc, argv);
    }
    if (command == "representative-ledger-manifest") {
      return representative_ledger_manifest(argc, argv);
    }
    if (command == "manifest-filter-quotient") {
      return filter_manifest_quotient_dimension(argc, argv);
    }
    if (command == "manifest-support-canonical-keys") {
      return canonicalize_manifest_supports(argc, argv);
    }
    if (command == "manifest-merge") return merge_manifests(argc, argv);
    if (command == "sparse-completion-gate") return sparse_completion_gate();
    if (command == "sparse-decode") return sparse_decode(argc, argv);
    if (command == "direct-canonical") {
      return direct_tensor_canonicalize(argc, argv);
    }
    if (command == "direct-canonical-labels") {
      return direct_tensor_canonicalize_labels(argc, argv);
    }
    if (command == "direct-tensor-gate") return direct_tensor_gate();
    if (command == "five-qubit-tensor-orbits") {
      return five_qubit_tensor_orbit_census(argc, argv);
    }
    if (command == "tensor-extension-census") {
      return tensor_extension_census(argc, argv);
    }
    if (command == "tensor-unmark-q7") {
      return tensor_unmark_q7(argc, argv);
    }
    if (command == "tensor-centroid") return tensor_centroid(argc, argv);
    if (command == "decomposed-canonical") {
      return decomposed_tensor_canonicalize(argc, argv);
    }
    if (command == "hereditary-punctures") {
      return hereditary_punctures(argc, argv);
    }
    if (command == "origin-orbits") return origin_orbits(argc, argv);
    usage();
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "utsp-native: " << error.what() << '\n';
    return 1;
  }
}
