#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr const char* kMagic = "L48_M08_NATIVE_TASK_V1";
constexpr int kPointCount = 128;
constexpr int kMaximumLiftDimension = 18;

struct Mask128 {
  std::uint64_t low = 0;
  std::uint64_t high = 0;

  bool operator==(const Mask128&) const = default;
};

struct Mask128Hash {
  std::size_t operator()(const Mask128& value) const noexcept {
    std::uint64_t mixed = value.low ^
        (value.high + 0x9e3779b97f4a7c15ULL + (value.low << 6) +
         (value.low >> 2));
    mixed ^= mixed >> 30;
    mixed *= 0xbf58476d1ce4e5b9ULL;
    mixed ^= mixed >> 27;
    mixed *= 0x94d049bb133111ebULL;
    mixed ^= mixed >> 31;
    return static_cast<std::size_t>(mixed);
  }
};

bool mask_less(const Mask128& left, const Mask128& right) {
  return left.high < right.high ||
         (left.high == right.high && left.low < right.low);
}

void set_point(Mask128& mask, int point) {
  if (point < 0 || point >= kPointCount) {
    throw std::runtime_error("Point label outside F_2^7");
  }
  if (point < 64) {
    mask.low |= 1ULL << point;
  } else {
    mask.high |= 1ULL << (point - 64);
  }
}

std::uint64_t degree_two_vector(int point) {
  std::array<int, 7> coordinates{};
  for (int variable = 0; variable < 7; ++variable) {
    coordinates[variable] = (point >> (6 - variable)) & 1;
  }
  std::uint64_t value = 1;
  std::uint64_t bit = 1;
  for (const int coordinate : coordinates) {
    bit <<= 1;
    if (coordinate != 0) value |= bit;
  }
  for (int left = 0; left < 7; ++left) {
    for (int right = left + 1; right < 7; ++right) {
      bit <<= 1;
      if (coordinates[left] != 0 && coordinates[right] != 0) value |= bit;
    }
  }
  return value;
}

template <typename Visitor>
void visit_points(Mask128 mask, Visitor visitor) {
  while (mask.low != 0) {
    const int bit = __builtin_ctzll(mask.low);
    visitor(bit);
    mask.low &= mask.low - 1;
  }
  while (mask.high != 0) {
    const int bit = __builtin_ctzll(mask.high);
    visitor(bit + 64);
    mask.high &= mask.high - 1;
  }
}

struct CombinationSolver {
  std::array<std::uint64_t, 64> rows{};
  std::array<std::uint64_t, 64> coefficients{};
  int column_count = 0;
  int rank = 0;

  explicit CombinationSolver(const std::vector<std::uint64_t>& columns) {
    if (columns.size() > 64) {
      throw std::runtime_error("Combination solver has more than 64 columns");
    }
    column_count = static_cast<int>(columns.size());
    for (int index = 0; index < column_count; ++index) {
      std::uint64_t value = columns[index];
      std::uint64_t combination = 1ULL << index;
      for (int pivot = 63; pivot >= 0; --pivot) {
        if (((value >> pivot) & 1ULL) != 0 && rows[pivot] != 0) {
          value ^= rows[pivot];
          combination ^= coefficients[pivot];
        }
      }
      if (value == 0) continue;
      const int pivot = 63 - __builtin_clzll(value);
      rows[pivot] = value;
      coefficients[pivot] = combination;
      ++rank;
    }
  }

  std::optional<std::uint64_t> solve(std::uint64_t target) const {
    std::uint64_t combination = 0;
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (((target >> pivot) & 1ULL) != 0 && rows[pivot] != 0) {
        target ^= rows[pivot];
        combination ^= coefficients[pivot];
      }
    }
    if (target != 0) return std::nullopt;
    return combination;
  }
};

struct AffineAction {
  std::uint32_t translation = 0;
  std::array<std::uint32_t, kMaximumLiftDimension> images{};
  std::uint8_t dimension = 0;

  bool operator==(const AffineAction& other) const {
    if (translation != other.translation || dimension != other.dimension) {
      return false;
    }
    for (int index = 0; index < dimension; ++index) {
      if (images[index] != other.images[index]) return false;
    }
    return true;
  }
};

struct AffineActionHash {
  std::size_t operator()(const AffineAction& action) const noexcept {
    std::uint64_t hash = action.translation + 0x9e3779b97f4a7c15ULL;
    hash ^= static_cast<std::uint64_t>(action.dimension) << 56;
    for (int index = 0; index < action.dimension; ++index) {
      hash ^= static_cast<std::uint64_t>(action.images[index]) +
              0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return static_cast<std::size_t>(hash);
  }
};

bool action_less(const AffineAction& left, const AffineAction& right) {
  if (left.translation != right.translation) {
    return left.translation < right.translation;
  }
  for (int index = 0; index < left.dimension; ++index) {
    if (left.images[index] != right.images[index]) {
      return left.images[index] < right.images[index];
    }
  }
  return false;
}

AffineAction identity_action(int dimension) {
  AffineAction result;
  result.dimension = static_cast<std::uint8_t>(dimension);
  for (int index = 0; index < dimension; ++index) {
    result.images[index] = 1U << index;
  }
  return result;
}

std::uint32_t linear_apply(const AffineAction& action, std::uint32_t value) {
  std::uint32_t result = 0;
  while (value != 0) {
    const int bit = __builtin_ctz(value);
    result ^= action.images[bit];
    value &= value - 1;
  }
  return result;
}

std::uint32_t apply(const AffineAction& action, std::uint32_t value) {
  return action.translation ^ linear_apply(action, value);
}

AffineAction compose(const AffineAction& first, const AffineAction& second) {
  if (first.dimension != second.dimension) {
    throw std::runtime_error("Affine action dimension mismatch");
  }
  AffineAction result;
  result.dimension = first.dimension;
  result.translation = apply(second, first.translation);
  for (int index = 0; index < result.dimension; ++index) {
    result.images[index] = linear_apply(second, first.images[index]);
  }
  return result;
}

AffineAction inverse(const AffineAction& action) {
  std::vector<std::uint64_t> columns;
  for (int index = 0; index < action.dimension; ++index) {
    columns.push_back(action.images[index]);
  }
  CombinationSolver solver(columns);
  if (solver.rank != action.dimension) {
    throw std::runtime_error("Cannot invert singular affine action");
  }
  AffineAction result;
  result.dimension = action.dimension;
  for (int index = 0; index < action.dimension; ++index) {
    const auto solution = solver.solve(1ULL << index);
    if (!solution.has_value()) {
      throw std::runtime_error("Failed to invert affine action");
    }
    result.images[index] = static_cast<std::uint32_t>(*solution);
  }
  result.translation = linear_apply(result, action.translation);
  return result;
}

struct Task {
  std::string parent_id;
  int weight = 0;
  int multiplicity = 0;
  int core_count = 0;
  int external_count = 0;
  int affine_count = 0;
  int lift_dimension = 0;
  int generator_count = 0;
  std::uint64_t stabilizer_order = 0;
  std::uint64_t expected_fiber_count = 0;
  std::uint64_t expected_raw_count = 0;
  std::uint64_t expected_excluded_count = 0;
  std::vector<int> core_points;
  std::vector<std::uint64_t> quadratic_columns;
  std::vector<int> external_points;
  std::vector<std::uint64_t> external_labels;
  std::vector<std::uint64_t> affine_rows;
  std::vector<std::uint64_t> lift_complement;
  std::vector<std::array<std::uint8_t, kPointCount>> generators;
};

template <typename Value>
std::vector<Value> read_vector(std::istream& input, int count) {
  std::vector<Value> result(count);
  for (int index = 0; index < count; ++index) {
    unsigned long long value = 0;
    input >> value;
    if (!input) throw std::runtime_error("Truncated native task input");
    result[index] = static_cast<Value>(value);
  }
  return result;
}

Task read_task(const fs::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Could not open native task input");
  std::string magic;
  input >> magic;
  if (magic != kMagic) throw std::runtime_error("Native task magic mismatch");
  Task task;
  input >> task.parent_id >> task.weight >> task.multiplicity >> task.core_count >>
      task.external_count >> task.affine_count >> task.lift_dimension >>
      task.generator_count >> task.stabilizer_order >> task.expected_fiber_count >>
      task.expected_raw_count >> task.expected_excluded_count;
  if (!input || task.weight + 2 * task.multiplicity != 48 ||
      task.core_count != task.weight || (task.affine_count != 7 && task.affine_count != 8) ||
      task.lift_dimension < 0 || task.lift_dimension > kMaximumLiftDimension ||
      task.generator_count <= 0 || task.stabilizer_order == 0) {
    throw std::runtime_error("Invalid native task header");
  }
  task.core_points = read_vector<int>(input, task.core_count);
  task.quadratic_columns =
      read_vector<std::uint64_t>(input, task.core_count);
  task.external_points = read_vector<int>(input, task.external_count);
  task.external_labels =
      read_vector<std::uint64_t>(input, task.external_count);
  task.affine_rows = read_vector<std::uint64_t>(input, task.affine_count);
  task.lift_complement =
      read_vector<std::uint64_t>(input, task.lift_dimension);
  task.generators.resize(task.generator_count);
  for (auto& generator : task.generators) {
    const auto values = read_vector<int>(input, kPointCount);
    std::array<bool, kPointCount> seen{};
    for (int point = 0; point < kPointCount; ++point) {
      if (values[point] < 0 || values[point] >= kPointCount ||
          seen[values[point]]) {
        throw std::runtime_error("A supplied generator is not a permutation");
      }
      seen[values[point]] = true;
      generator[point] = static_cast<std::uint8_t>(values[point]);
    }
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("Unexpected trailing task data");
  return task;
}

Mask128 fiber_mask(const std::vector<int>& points,
                   std::initializer_list<int> indices) {
  Mask128 result;
  for (const int index : indices) set_point(result, points[index]);
  return result;
}

std::vector<Mask128> generate_valid_fibers(const Task& task) {
  std::vector<Mask128> result;
  if (task.multiplicity == 0) {
    result.push_back({});
  } else if (task.multiplicity == 2) {
    for (int left = 0; left < task.external_count - 1; ++left) {
      for (int right = left + 1; right < task.external_count; ++right) {
        if (task.external_labels[left] == task.external_labels[right]) {
          result.push_back(fiber_mask(task.external_points, {left, right}));
        }
      }
    }
  } else if (task.multiplicity == 4) {
    std::map<std::uint64_t, std::vector<std::pair<int, int>>> pairs_by_xor;
    for (int left = 0; left < task.external_count - 1; ++left) {
      for (int right = left + 1; right < task.external_count; ++right) {
        pairs_by_xor[task.external_labels[left] ^ task.external_labels[right]]
            .emplace_back(left, right);
      }
    }
    for (const auto& [label, pairs] : pairs_by_xor) {
      (void)label;
      for (std::size_t first = 0; first + 1 < pairs.size(); ++first) {
        for (std::size_t second = first + 1; second < pairs.size(); ++second) {
          if (pairs[first].second < pairs[second].first) {
            result.push_back(fiber_mask(
                task.external_points,
                {pairs[first].first, pairs[first].second,
                 pairs[second].first, pairs[second].second}));
          }
        }
      }
    }
  } else {
    throw std::runtime_error("Unsupported complete-fibre multiplicity");
  }
  std::sort(result.begin(), result.end(), mask_less);
  result.erase(std::unique(result.begin(), result.end()), result.end());
  if (result.size() != task.expected_fiber_count) {
    throw std::runtime_error("Native valid-fibre count disagrees with theory gate");
  }
  return result;
}

Mask128 permute_fiber(
    Mask128 source,
    const std::array<std::uint8_t, kPointCount>& generator) {
  Mask128 result;
  visit_points(source, [&](int point) { set_point(result, generator[point]); });
  return result;
}

std::uint64_t permute_core_mask(
    std::uint64_t source,
    const std::vector<std::uint8_t>& permutation) {
  std::uint64_t result = 0;
  while (source != 0) {
    const int bit = __builtin_ctzll(source);
    result |= 1ULL << permutation[bit];
    source &= source - 1;
  }
  return result;
}

struct OutputRecord {
  std::array<std::uint64_t, 4> support{};
  std::uint64_t marked_orbit_size = 0;
};

static_assert(sizeof(OutputRecord) == 40);

void set_support_point(std::array<std::uint64_t, 4>& support, int point) {
  if (point < 0 || point >= 256) {
    throw std::runtime_error("Support point outside F_2^8");
  }
  support[point / 64] |= 1ULL << (point % 64);
}

int support_weight(const std::array<std::uint64_t, 4>& support) {
  int result = 0;
  for (const auto word : support) result += __builtin_popcountll(word);
  return result;
}

struct Arguments {
  fs::path input;
  fs::path output;
  fs::path summary;
  fs::path debug;
  bool count_only = false;
};

Arguments parse_arguments(int argc, char** argv) {
  Arguments result;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (option == "--count-only") {
      result.count_only = true;
      continue;
    }
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const fs::path value = argv[++index];
    if (option == "--input") {
      result.input = value;
    } else if (option == "--output") {
      result.output = value;
    } else if (option == "--summary") {
      result.summary = value;
    } else if (option == "--debug") {
      result.debug = value;
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (result.input.empty() || result.summary.empty() ||
      (!result.count_only && result.output.empty()) ||
      (result.count_only && !result.output.empty())) {
    throw std::runtime_error(
        "Usage: --input TASK --summary JSON [--count-only | --output BIN]");
  }
  return result;
}

class Kernel {
 public:
  Kernel(Task task, const Arguments& arguments)
      : task_(std::move(task)), arguments_(arguments),
        quadratic_solver_(task_.quadratic_columns) {
    std::vector<std::uint64_t> kernel_columns = task_.affine_rows;
    kernel_columns.insert(kernel_columns.end(), task_.lift_complement.begin(),
                          task_.lift_complement.end());
    kernel_solver_.emplace(kernel_columns);
    if (quadratic_solver_.rank + task_.affine_count + task_.lift_dimension !=
            task_.core_count ||
        kernel_solver_->rank != task_.affine_count + task_.lift_dimension) {
      throw std::runtime_error("Task rank-nullity data are inconsistent");
    }
    for (int point = 0; point < kPointCount; ++point) core_index_[point] = -1;
    for (int index = 0; index < task_.core_count; ++index) {
      const int point = task_.core_points[index];
      if (point < 0 || point >= kPointCount || core_index_[point] != -1) {
        throw std::runtime_error("Core support is not projective");
      }
      core_index_[point] = index;
    }
    for (int point = 0; point < kPointCount; ++point) {
      quadratic_by_point_[point] = degree_two_vector(point);
    }
    build_generator_actions();
  }

  void run() {
    valid_fibers_ = generate_valid_fibers(task_);
    fiber_index_.reserve(valid_fibers_.size() * 2);
    for (std::uint32_t index = 0; index < valid_fibers_.size(); ++index) {
      fiber_index_.emplace(valid_fibers_[index], index);
    }
    particulars_.resize(valid_fibers_.size());
    for (std::size_t index = 0; index < valid_fibers_.size(); ++index) {
      std::uint64_t syndrome = 0;
      visit_points(valid_fibers_[index], [&](int point) {
        syndrome ^= quadratic_by_point_[point];
      });
      const auto particular = quadratic_solver_.solve(syndrome);
      if (!particular.has_value()) {
        throw std::runtime_error("Valid fibre has no particular singleton lift");
      }
      particulars_[index] = *particular;
      std::uint64_t evaluated = 0;
      std::uint64_t coefficients = *particular;
      while (coefficients != 0) {
        const int bit = __builtin_ctzll(coefficients);
        evaluated ^= task_.quadratic_columns[bit];
        coefficients &= coefficients - 1;
      }
      if (evaluated != syndrome) {
        throw std::runtime_error("Particular singleton lift failed evaluation");
      }
    }

    if (!arguments_.count_only) {
      if (fs::exists(arguments_.output)) {
        throw std::runtime_error("Refusing to overwrite native candidate output");
      }
      output_temporary_ = arguments_.output.string() + ".tmp";
      output_.open(output_temporary_, std::ios::binary | std::ios::trunc);
      if (!output_) throw std::runtime_error("Could not open candidate output");
    }
    if (!arguments_.debug.empty()) {
      if (fs::exists(arguments_.debug)) {
        throw std::runtime_error("Refusing to overwrite native debug output");
      }
      debug_temporary_ = arguments_.debug.string() + ".tmp";
      debug_.open(debug_temporary_, std::ios::trunc);
      if (!debug_) throw std::runtime_error("Could not open native debug output");
      debug_ << "seed_low\tseed_high\tfiber_orbit_size\tinduced_generator_count"
                "\tlift_orbit_count\tseed_particular\tlift_orbits"
                "\tinduced_actions\n";
    }

    std::vector<std::uint8_t> completed(valid_fibers_.size(), 0);
    positions_.assign(valid_fibers_.size(), -1);
    for (std::uint32_t seed = 0; seed < valid_fibers_.size(); ++seed) {
      if (completed[seed]) continue;
      process_fiber_orbit(seed, completed);
    }
    if (raw_accounting_ != task_.expected_raw_count ||
        excluded_count_ != task_.expected_excluded_count) {
      throw std::runtime_error("Native orbit accounting disagrees with theory gate");
    }
    if (!arguments_.count_only) {
      output_.flush();
      if (!output_) throw std::runtime_error("Candidate output write failed");
      output_.close();
      if (fs::file_size(output_temporary_) !=
          direction_marked_orbit_count_ * sizeof(OutputRecord)) {
        throw std::runtime_error("Candidate output has the wrong byte count");
      }
      fs::rename(output_temporary_, arguments_.output);
    }
    if (!arguments_.debug.empty()) {
      debug_.flush();
      if (!debug_) throw std::runtime_error("Native debug output write failed");
      debug_.close();
      fs::rename(debug_temporary_, arguments_.debug);
    }
  }

  void write_summary(double elapsed_seconds) const {
    if (fs::exists(arguments_.summary)) {
      throw std::runtime_error("Refusing to overwrite native summary");
    }
    const fs::path temporary = arguments_.summary.string() + ".tmp";
    std::ofstream output(temporary, std::ios::trunc);
    output << "{\n"
           << "  \"status\": \"complete_length48_m08_native_marked_orbits_v1\",\n"
           << "  \"parent_id\": \"" << task_.parent_id << "\",\n"
           << "  \"parent_weight\": " << task_.weight << ",\n"
           << "  \"full_fiber_multiplicity\": " << task_.multiplicity << ",\n"
           << "  \"compatible_fiber_set_count\": " << valid_fibers_.size()
           << ",\n"
           << "  \"fiber_orbit_count\": " << fiber_orbit_count_ << ",\n"
           << "  \"lift_quotient_dimension\": " << task_.lift_dimension << ",\n"
           << "  \"core_stabilizer_order\": " << task_.stabilizer_order << ",\n"
           << "  \"deduplicated_generator_count\": " << task_.generator_count
           << ",\n"
           << "  \"direction_marked_orbit_count\": "
           << direction_marked_orbit_count_ << ",\n"
           << "  \"raw_marked_pair_count\": " << raw_accounting_ << ",\n"
           << "  \"excluded_rank_deficient_pair_count\": " << excluded_count_
           << ",\n"
           << "  \"schreier_edge_count\": " << schreier_edge_count_ << ",\n"
           << "  \"maximum_fiber_orbit_size\": " << maximum_fiber_orbit_size_
           << ",\n"
           << "  \"candidate_binary_bytes\": "
           << (arguments_.count_only
                   ? 0
                   : direction_marked_orbit_count_ * sizeof(OutputRecord))
           << ",\n"
           << "  \"count_only\": "
           << (arguments_.count_only ? "true" : "false") << ",\n"
           << "  \"elapsed_seconds\": " << std::setprecision(17)
           << elapsed_seconds << ",\n"
           << "  \"all_checks_pass\": true\n"
           << "}\n";
    if (!output) throw std::runtime_error("Could not write native summary");
    output.close();
    fs::rename(temporary, arguments_.summary);
  }

  void cleanup() const {
    if (!output_temporary_.empty()) {
      std::error_code ignored;
      fs::remove(output_temporary_, ignored);
    }
    if (!debug_temporary_.empty()) {
      std::error_code ignored;
      fs::remove(debug_temporary_, ignored);
    }
  }

 private:
  void build_generator_actions() {
    const auto identity = identity_action(task_.lift_dimension);
    (void)identity;
    for (const auto& generator : task_.generators) {
      std::vector<std::uint8_t> core_permutation(task_.core_count);
      for (int index = 0; index < task_.core_count; ++index) {
        const int image = generator[task_.core_points[index]];
        if (core_index_[image] < 0) {
          throw std::runtime_error("A generator does not preserve the core");
        }
        core_permutation[index] =
            static_cast<std::uint8_t>(core_index_[image]);
      }
      AffineAction linear;
      linear.dimension = static_cast<std::uint8_t>(task_.lift_dimension);
      for (int index = 0; index < task_.lift_dimension; ++index) {
        const auto transformed =
            permute_core_mask(task_.lift_complement[index], core_permutation);
        const auto coefficients = kernel_solver_->solve(transformed);
        if (!coefficients.has_value()) {
          throw std::runtime_error("Generator left the lift kernel");
        }
        linear.images[index] = static_cast<std::uint32_t>(
            *coefficients >> task_.affine_count);
      }
      if (inverse(linear).dimension != task_.lift_dimension) {
        throw std::runtime_error("Singular generator action");
      }
      core_permutations_.push_back(std::move(core_permutation));
      linear_actions_.push_back(linear);
    }
  }

  AffineAction transition(std::uint32_t source, std::uint32_t target,
                          int generator_index) const {
    const auto transformed = permute_core_mask(
        particulars_[source], core_permutations_[generator_index]);
    const auto coefficients =
        kernel_solver_->solve(transformed ^ particulars_[target]);
    if (!coefficients.has_value()) {
      throw std::runtime_error("Chart transition left the lift kernel");
    }
    AffineAction result = linear_actions_[generator_index];
    result.translation = static_cast<std::uint32_t>(
        *coefficients >> task_.affine_count);
    return result;
  }

  std::uint32_t image_index(std::uint32_t source, int generator_index) const {
    const auto image =
        permute_fiber(valid_fibers_[source], task_.generators[generator_index]);
    const auto iterator = fiber_index_.find(image);
    if (iterator == fiber_index_.end()) {
      throw std::runtime_error("Core automorphism left valid-fibre domain");
    }
    return iterator->second;
  }

  std::vector<std::pair<std::uint32_t, std::uint32_t>> lift_orbits(
      const std::unordered_set<AffineAction, AffineActionHash>& generators) const {
    const auto identity = identity_action(task_.lift_dimension);
    std::unordered_set<AffineAction, AffineActionHash> closed;
    for (const auto& generator : generators) {
      if (!(generator == identity)) {
        closed.insert(generator);
        closed.insert(inverse(generator));
      }
    }
    std::vector<AffineAction> actions(closed.begin(), closed.end());
    std::sort(actions.begin(), actions.end(), action_less);
    const std::uint32_t size = 1U << task_.lift_dimension;
    if (actions.empty()) {
      std::vector<std::pair<std::uint32_t, std::uint32_t>> result;
      result.reserve(size);
      for (std::uint32_t value = 0; value < size; ++value) {
        result.emplace_back(value, 1);
      }
      return result;
    }
    std::vector<std::uint8_t> seen(size, 0);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> result;
    std::vector<std::uint32_t> queue;
    for (std::uint32_t seed = 0; seed < size; ++seed) {
      if (seen[seed]) continue;
      seen[seed] = 1;
      queue.clear();
      queue.push_back(seed);
      for (std::size_t offset = 0; offset < queue.size(); ++offset) {
        for (const auto& action : actions) {
          const auto image = apply(action, queue[offset]);
          if (!seen[image]) {
            seen[image] = 1;
            queue.push_back(image);
          }
        }
      }
      result.emplace_back(seed, static_cast<std::uint32_t>(queue.size()));
    }
    return result;
  }

  OutputRecord output_record(std::uint32_t fiber_index,
                             std::uint32_t lift_coordinates,
                             std::uint64_t orbit_size) const {
    std::uint64_t lift = particulars_[fiber_index];
    std::uint32_t remaining = lift_coordinates;
    while (remaining != 0) {
      const int bit = __builtin_ctz(remaining);
      lift ^= task_.lift_complement[bit];
      remaining &= remaining - 1;
    }
    OutputRecord record;
    for (int index = 0; index < task_.core_count; ++index) {
      const int point = 2 * task_.core_points[index] +
                        static_cast<int>((lift >> index) & 1ULL);
      set_support_point(record.support, point);
    }
    visit_points(valid_fibers_[fiber_index], [&](int point) {
      set_support_point(record.support, 2 * point);
      set_support_point(record.support, 2 * point + 1);
    });
    if (support_weight(record.support) != 48) {
      throw std::runtime_error("Reconstructed candidate has wrong weight");
    }
    record.marked_orbit_size = orbit_size;
    return record;
  }

  void process_fiber_orbit(std::uint32_t seed,
                           std::vector<std::uint8_t>& completed) {
    const auto identity = identity_action(task_.lift_dimension);
    std::vector<std::uint32_t> orbit;
    std::vector<AffineAction> transporters;
    std::vector<AffineAction> inverse_transporters;
    orbit.push_back(seed);
    transporters.push_back(identity);
    inverse_transporters.push_back(identity);
    positions_[seed] = 0;
    std::unordered_set<AffineAction, AffineActionHash> stabilizer_actions;

    for (std::size_t offset = 0; offset < orbit.size(); ++offset) {
      const auto current = orbit[offset];
      for (int generator = 0; generator < task_.generator_count; ++generator) {
        ++schreier_edge_count_;
        const auto image = image_index(current, generator);
        const auto chart = transition(current, image, generator);
        const auto candidate = compose(transporters[offset], chart);
        if (positions_[image] < 0) {
          positions_[image] = static_cast<std::int32_t>(orbit.size());
          orbit.push_back(image);
          transporters.push_back(candidate);
          inverse_transporters.push_back(inverse(candidate));
        } else {
          const auto schreier = compose(
              candidate, inverse_transporters[positions_[image]]);
          if (!(schreier == identity)) stabilizer_actions.insert(schreier);
        }
      }
    }
    if (task_.stabilizer_order % orbit.size() != 0) {
      throw std::runtime_error("Fibre orbit does not divide stabilizer order");
    }
    for (const auto index : orbit) {
      if (completed[index]) throw std::runtime_error("Fibre orbits overlap");
      completed[index] = 1;
      positions_[index] = -1;
    }
    ++fiber_orbit_count_;
    maximum_fiber_orbit_size_ =
        std::max<std::uint64_t>(maximum_fiber_orbit_size_, orbit.size());
    const auto computed_lift_orbits = lift_orbits(stabilizer_actions);
    if (debug_) {
      debug_ << valid_fibers_[seed].low << '\t' << valid_fibers_[seed].high
             << '\t' << orbit.size() << '\t' << stabilizer_actions.size()
             << '\t' << computed_lift_orbits.size() << '\t'
             << particulars_[seed] << '\t';
      for (std::size_t index = 0; index < computed_lift_orbits.size(); ++index) {
        if (index != 0) debug_ << ',';
        debug_ << computed_lift_orbits[index].first << ':'
               << computed_lift_orbits[index].second;
      }
      debug_ << '\t';
      std::vector<AffineAction> ordered_stabilizer_actions(
          stabilizer_actions.begin(), stabilizer_actions.end());
      std::sort(ordered_stabilizer_actions.begin(),
                ordered_stabilizer_actions.end(), action_less);
      for (std::size_t index = 0; index < ordered_stabilizer_actions.size();
           ++index) {
        if (index != 0) debug_ << '|';
        const auto& action = ordered_stabilizer_actions[index];
        debug_ << action.translation << ';';
        for (int column = 0; column < action.dimension; ++column) {
          if (column != 0) debug_ << ',';
          debug_ << action.images[column];
        }
      }
      debug_ << '\n';
    }
    for (const auto [lift, lift_orbit_size] : computed_lift_orbits) {
      const std::uint64_t marked_size = orbit.size() * lift_orbit_size;
      if (task_.multiplicity == 0 && lift == 0) {
        excluded_count_ += marked_size;
        continue;
      }
      raw_accounting_ += marked_size;
      ++direction_marked_orbit_count_;
      if (!arguments_.count_only) {
        const auto record = output_record(seed, lift, marked_size);
        output_.write(reinterpret_cast<const char*>(&record), sizeof(record));
        if (!output_) throw std::runtime_error("Candidate output write failed");
      }
    }
  }

  Task task_;
  Arguments arguments_;
  CombinationSolver quadratic_solver_;
  std::optional<CombinationSolver> kernel_solver_;
  std::array<int, kPointCount> core_index_{};
  std::array<std::uint64_t, kPointCount> quadratic_by_point_{};
  std::vector<std::vector<std::uint8_t>> core_permutations_;
  std::vector<AffineAction> linear_actions_;
  std::vector<Mask128> valid_fibers_;
  std::unordered_map<Mask128, std::uint32_t, Mask128Hash> fiber_index_;
  std::vector<std::uint64_t> particulars_;
  std::vector<std::int32_t> positions_;
  std::ofstream output_;
  std::ofstream debug_;
  fs::path output_temporary_;
  fs::path debug_temporary_;
  std::uint64_t fiber_orbit_count_ = 0;
  std::uint64_t direction_marked_orbit_count_ = 0;
  std::uint64_t raw_accounting_ = 0;
  std::uint64_t excluded_count_ = 0;
  std::uint64_t schreier_edge_count_ = 0;
  std::uint64_t maximum_fiber_orbit_size_ = 0;
};

}  // namespace

int main(int argc, char** argv) {
  std::optional<Kernel> kernel;
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto task = read_task(arguments.input);
    const auto started = std::chrono::steady_clock::now();
    kernel.emplace(task, arguments);
    kernel->run();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
    kernel->write_summary(elapsed);
    std::cout << "parent_id=" << task.parent_id << "\n"
              << "elapsed_seconds=" << std::setprecision(17) << elapsed << "\n"
              << "all_checks_pass=true\n";
    return 0;
  } catch (const std::exception& error) {
    if (kernel.has_value()) kernel->cleanup();
    std::cerr << "m08_marked_orbit_kernel: " << error.what() << '\n';
    return 1;
  }
}
