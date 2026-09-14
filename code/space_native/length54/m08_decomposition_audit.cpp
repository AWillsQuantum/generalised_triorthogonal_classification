#include <algorithm>
#include <array>
#include <bit>
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

constexpr int kTargetLength = 54;

struct RepresentativeRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 4> mask{};
  std::uint64_t member_count = 0;
};

struct AuditRecord {
  std::array<std::uint64_t, 4> mask{};
  std::uint64_t member_count = 0;
  std::uint32_t matroid_component_count = 0;
  std::uint32_t stabilizer_dimension = 0;
  std::uint32_t indecomposable = 0;
  std::uint32_t reserved = 0;
};

static_assert(sizeof(RepresentativeRecord) == 72);
static_assert(sizeof(AuditRecord) == 56);

struct Arguments {
  fs::path input;
  fs::path output;
  std::uint64_t expected_records = 0;
};

struct UnionFind {
  std::array<int, kTargetLength> parent{};
  std::array<int, kTargetLength> rank{};

  UnionFind() {
    for (int index = 0; index < kTargetLength; ++index) parent[index] = index;
  }
  int find(int value) {
    if (parent[value] != value) parent[value] = find(parent[value]);
    return parent[value];
  }
  void unite(int left, int right) {
    left = find(left);
    right = find(right);
    if (left == right) return;
    if (rank[left] < rank[right]) std::swap(left, right);
    parent[right] = left;
    if (rank[left] == rank[right]) ++rank[left];
  }
};

std::vector<int> support_points(const std::array<std::uint64_t, 4>& mask) {
  std::vector<int> points;
  points.reserve(kTargetLength);
  for (int block = 0; block < 4; ++block) {
    auto word = mask[block];
    while (word != 0) {
      const int bit = std::countr_zero(word);
      points.push_back(64 * block + bit);
      word &= word - 1;
    }
  }
  return points;
}

int binary_rank(std::vector<std::uint64_t> values) {
  int rank = 0;
  for (int bit = 63; bit >= 0; --bit) {
    auto pivot = std::find_if(
        values.begin() + rank, values.end(),
        [bit](std::uint64_t value) { return ((value >> bit) & 1U) != 0; });
    if (pivot == values.end()) continue;
    std::iter_swap(values.begin() + rank, pivot);
    for (std::size_t row = 0; row < values.size(); ++row) {
      if (static_cast<int>(row) != rank &&
          ((values[row] >> bit) & 1U) != 0) {
        values[row] ^= values[rank];
      }
    }
    ++rank;
  }
  return rank;
}

int matroid_component_count(const std::vector<int>& points) {
  if (points.size() != kTargetLength) {
    throw std::runtime_error("Support weight is not 54");
  }
  std::vector<std::uint64_t> columns;
  columns.reserve(kTargetLength);
  for (const int point : points) {
    columns.push_back(static_cast<std::uint64_t>(1U | (point << 1)));
  }
  const int target_rank = binary_rank(columns);
  std::vector<int> basis_indices;
  std::vector<std::uint64_t> basis_columns;
  for (int index = 0; index < kTargetLength; ++index) {
    auto candidate = basis_columns;
    candidate.push_back(columns[index]);
    if (binary_rank(candidate) > static_cast<int>(basis_columns.size())) {
      basis_indices.push_back(index);
      basis_columns.push_back(columns[index]);
    }
    if (static_cast<int>(basis_columns.size()) == target_rank) break;
  }
  if (target_rank != 9 || basis_columns.size() != 9) {
    throw std::runtime_error("Representative does not have affine rank eight");
  }

  std::array<std::uint64_t, 64> echelon{};
  std::array<std::uint64_t, 64> coefficients{};
  std::array<bool, 64> occupied{};
  for (int basis_position = 0; basis_position < 9; ++basis_position) {
    auto value = basis_columns[basis_position];
    std::uint64_t coeff = 1ULL << basis_position;
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (occupied[pivot] && ((value >> pivot) & 1U) != 0) {
        value ^= echelon[pivot];
        coeff ^= coefficients[pivot];
      }
    }
    if (value == 0) throw std::runtime_error("Chosen basis is dependent");
    const int pivot = 63 - std::countl_zero(value);
    occupied[pivot] = true;
    echelon[pivot] = value;
    coefficients[pivot] = coeff;
  }

  std::array<bool, kTargetLength> is_basis{};
  for (const int index : basis_indices) is_basis[index] = true;
  UnionFind components;
  for (int index = 0; index < kTargetLength; ++index) {
    if (is_basis[index]) continue;
    auto value = columns[index];
    std::uint64_t coeff = 0;
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (occupied[pivot] && ((value >> pivot) & 1U) != 0) {
        value ^= echelon[pivot];
        coeff ^= coefficients[pivot];
      }
    }
    if (value != 0) throw std::runtime_error("Column lies outside basis span");
    for (int basis_position = 0; basis_position < 9; ++basis_position) {
      if (((coeff >> basis_position) & 1U) != 0) {
        components.unite(index, basis_indices[basis_position]);
      }
    }
  }
  std::array<bool, kTargetLength> roots{};
  for (int index = 0; index < kTargetLength; ++index) {
    roots[components.find(index)] = true;
  }
  return static_cast<int>(std::count(roots.begin(), roots.end(), true));
}

int multiplicative_stabilizer_dimension(const std::vector<int>& points) {
  // A multiplier is a length-54 vector lambda for which lambda*x_j is affine
  // on the support for every affine coordinate row x_j (including x_0=1).
  // For m=8 this is a 486-equation binary system in 54 unknowns.
  std::array<std::uint64_t, kTargetLength> coordinate_rows{};
  coordinate_rows[0] = (1ULL << kTargetLength) - 1;
  for (int bit = 0; bit < 8; ++bit) {
    for (int column = 0; column < kTargetLength; ++column) {
      if (((points[column] >> bit) & 1U) != 0) {
        coordinate_rows[bit + 1] |= 1ULL << column;
      }
    }
  }

  // Compute a basis for the dual of the 9-dimensional affine row space.
  std::vector<std::uint64_t> row_basis(
      coordinate_rows.begin(), coordinate_rows.begin() + 9);
  std::vector<int> pivots;
  int pivot_row = 0;
  for (int column = 0; column < kTargetLength && pivot_row < 9; ++column) {
    int candidate = -1;
    for (int row = pivot_row; row < 9; ++row) {
      if (((row_basis[row] >> column) & 1U) != 0) {
        candidate = row;
        break;
      }
    }
    if (candidate < 0) continue;
    std::swap(row_basis[pivot_row], row_basis[candidate]);
    for (int row = 0; row < 9; ++row) {
      if (row != pivot_row && ((row_basis[row] >> column) & 1U) != 0) {
        row_basis[row] ^= row_basis[pivot_row];
      }
    }
    pivots.push_back(column);
    ++pivot_row;
  }
  if (pivot_row != 9) throw std::runtime_error("Affine row space has wrong rank");
  std::array<bool, kTargetLength> is_pivot{};
  for (const int pivot : pivots) is_pivot[pivot] = true;
  std::vector<std::uint64_t> dual;
  dual.reserve(kTargetLength - 9);
  for (int free_column = 0; free_column < kTargetLength; ++free_column) {
    if (is_pivot[free_column]) continue;
    std::uint64_t vector = 1ULL << free_column;
    for (int row = 0; row < 9; ++row) {
      if (((row_basis[row] >> free_column) & 1U) != 0) {
        vector |= 1ULL << pivots[row];
      }
    }
    dual.push_back(vector);
  }

  std::vector<std::uint64_t> equations;
  equations.reserve(9 * dual.size());
  for (int row = 0; row < 9; ++row) {
    for (const auto dual_row : dual) {
      equations.push_back(coordinate_rows[row] & dual_row);
    }
  }
  return kTargetLength - binary_rank(std::move(equations));
}

Arguments parse_arguments(int argc, char** argv) {
  Arguments arguments;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (index + 1 >= argc) throw std::runtime_error("Missing option value");
    const std::string value = argv[++index];
    if (option == "--input") {
      arguments.input = value;
    } else if (option == "--output") {
      arguments.output = value;
    } else if (option == "--expected-records") {
      arguments.expected_records = std::stoull(value);
    } else {
      throw std::runtime_error("Unknown option: " + option);
    }
  }
  if (arguments.input.empty() || arguments.output.empty() ||
      arguments.expected_records == 0) {
    throw std::runtime_error(
        "Required: --input PATH --output PATH --expected-records N");
  }
  return arguments;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const auto arguments = parse_arguments(argc, argv);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(arguments.input) ||
        fs::file_size(arguments.input) !=
            arguments.expected_records * sizeof(RepresentativeRecord)) {
      throw std::runtime_error("Representative ledger has the wrong byte length");
    }
    if (fs::exists(arguments.output)) {
      throw std::runtime_error("Refusing to overwrite decomposition audit");
    }
    const auto temporary = arguments.output.string() + ".tmp";
    std::ifstream input(arguments.input, std::ios::binary);
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    std::uint64_t indecomposable_count = 0;
    std::uint64_t decomposable_count = 0;
    std::uint64_t method_disagreement_count = 0;
    for (std::uint64_t index = 0; index < arguments.expected_records; ++index) {
      RepresentativeRecord source;
      input.read(reinterpret_cast<char*>(&source), sizeof(source));
      if (!input) throw std::runtime_error("Representative ledger is truncated");
      const auto points = support_points(source.mask);
      const int matroid_components = matroid_component_count(points);
      const int stabilizer_dimension = multiplicative_stabilizer_dimension(points);
      const bool matroid_indecomposable = matroid_components == 1;
      const bool stabilizer_indecomposable = stabilizer_dimension == 1;
      if (matroid_indecomposable != stabilizer_indecomposable) {
        ++method_disagreement_count;
      }
      const bool indecomposable =
          matroid_indecomposable && stabilizer_indecomposable;
      indecomposable_count += indecomposable;
      decomposable_count += !indecomposable;
      const AuditRecord audit{
          source.mask,
          source.member_count,
          static_cast<std::uint32_t>(matroid_components),
          static_cast<std::uint32_t>(stabilizer_dimension),
          static_cast<std::uint32_t>(indecomposable),
          0,
      };
      output.write(reinterpret_cast<const char*>(&audit), sizeof(audit));
      if (!output) throw std::runtime_error("Could not write decomposition audit");
    }
    output.close();
    if (!output || method_disagreement_count != 0 ||
        fs::file_size(temporary) != arguments.expected_records * sizeof(AuditRecord)) {
      fs::remove(temporary);
      throw std::runtime_error("Decomposition audit validation failed");
    }
    fs::rename(temporary, arguments.output);
    const auto finished = std::chrono::steady_clock::now();
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m08_native_decomposition_audit\",\n"
              << "  \"representative_count\": " << arguments.expected_records << ",\n"
              << "  \"indecomposable_count\": " << indecomposable_count << ",\n"
              << "  \"decomposable_count\": " << decomposable_count << ",\n"
              << "  \"method_disagreement_count\": " << method_disagreement_count << ",\n"
              << "  \"elapsed_seconds\": "
              << std::chrono::duration<double>(finished - started).count() << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
