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

constexpr int kAffineDimension = 13;
constexpr int kCodeDimension = kAffineDimension + 1;
constexpr int kSupportSize = 54;

struct RepresentativeRecord {
  std::array<std::uint64_t, 4> signature{};
  std::array<std::uint64_t, 128> mask{};
  std::uint64_t member_count = 0;
};

struct AuditRecord {
  std::array<std::uint64_t, 128> mask{};
  std::uint64_t member_count = 0;
  std::uint32_t matroid_component_count = 0;
  std::uint32_t stabilizer_dimension = 0;
  std::uint32_t indecomposable = 0;
  std::uint32_t reserved = 0;
};

static_assert(sizeof(RepresentativeRecord) == 1064);
static_assert(sizeof(AuditRecord) == 1048);

struct UnionFind {
  std::array<int, kSupportSize> parent{};
  std::array<int, kSupportSize> rank{};
  UnionFind() {
    for (int index = 0; index < kSupportSize; ++index) parent[index] = index;
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

std::vector<int> support_points(
    const std::array<std::uint64_t, 128>& mask) {
  std::vector<int> points;
  points.reserve(kSupportSize);
  for (int block = 0; block < 128; ++block) {
    std::uint64_t word = mask[block];
    while (word != 0) {
      const int bit = __builtin_ctzll(word);
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
      if (static_cast<int>(row) != rank && ((values[row] >> bit) & 1U)) {
        values[row] ^= values[rank];
      }
    }
    ++rank;
  }
  return rank;
}

int matroid_component_count(const std::vector<int>& points) {
  if (points.size() != kSupportSize) {
    throw std::runtime_error("Support weight is not 54");
  }
  std::vector<std::uint64_t> columns;
  columns.reserve(kSupportSize);
  for (const int point : points) {
    columns.push_back(1ULL | (static_cast<std::uint64_t>(point) << 1));
  }
  const int target_rank = binary_rank(columns);
  std::vector<int> basis_indices;
  std::vector<std::uint64_t> basis_columns;
  for (int index = 0; index < kSupportSize; ++index) {
    auto candidate = basis_columns;
    candidate.push_back(columns[index]);
    if (binary_rank(candidate) > static_cast<int>(basis_columns.size())) {
      basis_indices.push_back(index);
      basis_columns.push_back(columns[index]);
    }
    if (basis_columns.size() == kCodeDimension) break;
  }
  if (target_rank != kCodeDimension || basis_columns.size() != kCodeDimension) {
    throw std::runtime_error("Representative has the wrong affine rank");
  }
  std::array<std::uint64_t, 64> echelon{};
  std::array<std::uint64_t, 64> coefficients{};
  std::array<bool, 64> occupied{};
  for (int basis_position = 0; basis_position < kCodeDimension;
       ++basis_position) {
    auto value = basis_columns[basis_position];
    std::uint64_t coefficient = 1ULL << basis_position;
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (occupied[pivot] && ((value >> pivot) & 1U)) {
        value ^= echelon[pivot];
        coefficient ^= coefficients[pivot];
      }
    }
    if (!value) throw std::runtime_error("Chosen basis is dependent");
    const int pivot = 63 - __builtin_clzll(value);
    occupied[pivot] = true;
    echelon[pivot] = value;
    coefficients[pivot] = coefficient;
  }
  std::array<bool, kSupportSize> is_basis{};
  for (const int index : basis_indices) is_basis[index] = true;
  UnionFind components;
  for (int index = 0; index < kSupportSize; ++index) {
    if (is_basis[index]) continue;
    auto value = columns[index];
    std::uint64_t coefficient = 0;
    for (int pivot = 63; pivot >= 0; --pivot) {
      if (occupied[pivot] && ((value >> pivot) & 1U)) {
        value ^= echelon[pivot];
        coefficient ^= coefficients[pivot];
      }
    }
    if (value) throw std::runtime_error("Column lies outside the basis span");
    for (int basis_position = 0; basis_position < kCodeDimension;
         ++basis_position) {
      if ((coefficient >> basis_position) & 1U) {
        components.unite(index, basis_indices[basis_position]);
      }
    }
  }
  std::array<bool, kSupportSize> roots{};
  for (int index = 0; index < kSupportSize; ++index) {
    roots[components.find(index)] = true;
  }
  return static_cast<int>(std::count(roots.begin(), roots.end(), true));
}

int multiplicative_stabilizer_dimension(const std::vector<int>& points) {
  std::array<std::uint64_t, kCodeDimension> coordinate_rows{};
  coordinate_rows[0] = (1ULL << kSupportSize) - 1;
  for (int bit = 0; bit < kAffineDimension; ++bit) {
    for (int column = 0; column < kSupportSize; ++column) {
      if ((points[column] >> bit) & 1U) {
        coordinate_rows[bit + 1] |= 1ULL << column;
      }
    }
  }
  std::vector<std::uint64_t> row_basis(coordinate_rows.begin(),
                                       coordinate_rows.end());
  std::vector<int> pivots;
  int pivot_row = 0;
  for (int column = 0; column < kSupportSize && pivot_row < kCodeDimension;
       ++column) {
    int candidate = -1;
    for (int row = pivot_row; row < kCodeDimension; ++row) {
      if ((row_basis[row] >> column) & 1U) {
        candidate = row;
        break;
      }
    }
    if (candidate < 0) continue;
    std::swap(row_basis[pivot_row], row_basis[candidate]);
    for (int row = 0; row < kCodeDimension; ++row) {
      if (row != pivot_row && ((row_basis[row] >> column) & 1U)) {
        row_basis[row] ^= row_basis[pivot_row];
      }
    }
    pivots.push_back(column);
    ++pivot_row;
  }
  if (pivot_row != kCodeDimension) {
    throw std::runtime_error("Affine row space has wrong rank");
  }
  std::array<bool, kSupportSize> is_pivot{};
  for (const int pivot : pivots) is_pivot[pivot] = true;
  std::vector<std::uint64_t> dual;
  for (int free_column = 0; free_column < kSupportSize; ++free_column) {
    if (is_pivot[free_column]) continue;
    std::uint64_t vector = 1ULL << free_column;
    for (int row = 0; row < kCodeDimension; ++row) {
      if ((row_basis[row] >> free_column) & 1U) {
        vector |= 1ULL << pivots[row];
      }
    }
    dual.push_back(vector);
  }
  std::vector<std::uint64_t> equations;
  equations.reserve(kCodeDimension * dual.size());
  for (int row = 0; row < kCodeDimension; ++row) {
    for (const auto dual_row : dual) {
      equations.push_back(coordinate_rows[row] & dual_row);
    }
  }
  return kSupportSize - binary_rank(std::move(equations));
}

}  // namespace

int main(int argc, char** argv) {
  fs::path temporary;
  try {
    if (argc != 4) {
      throw std::runtime_error("Required: INPUT OUTPUT EXPECTED_RECORDS");
    }
    const fs::path input_path = argv[1];
    const fs::path output_path = argv[2];
    const std::uint64_t expected_records = std::stoull(argv[3]);
    const auto started = std::chrono::steady_clock::now();
    if (!fs::is_regular_file(input_path) ||
        fs::file_size(input_path) !=
            expected_records * sizeof(RepresentativeRecord) ||
        fs::exists(output_path)) {
      throw std::runtime_error("Representative input or audit output is invalid");
    }
    temporary = output_path.string() + ".tmp";
    std::ifstream input(input_path, std::ios::binary);
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    std::uint64_t indecomposable_count = 0;
    std::uint64_t decomposable_count = 0;
    std::uint64_t disagreement_count = 0;
    for (std::uint64_t index = 0; index < expected_records; ++index) {
      RepresentativeRecord source;
      input.read(reinterpret_cast<char*>(&source), sizeof(source));
      if (!input) throw std::runtime_error("Representative ledger is truncated");
      const auto support = support_points(source.mask);
      const int components = matroid_component_count(support);
      const int stabilizer = multiplicative_stabilizer_dimension(support);
      const bool matroid_indec = components == 1;
      const bool stabilizer_indec = stabilizer == 1;
      disagreement_count += matroid_indec != stabilizer_indec;
      const bool indecomposable = matroid_indec && stabilizer_indec;
      indecomposable_count += indecomposable;
      decomposable_count += !indecomposable;
      const AuditRecord audit{
          source.mask,
          source.member_count,
          static_cast<std::uint32_t>(components),
          static_cast<std::uint32_t>(stabilizer),
          static_cast<std::uint32_t>(indecomposable),
          0};
      output.write(reinterpret_cast<const char*>(&audit), sizeof(audit));
    }
    output.close();
    if (!output || disagreement_count ||
        fs::file_size(temporary) != expected_records * sizeof(AuditRecord)) {
      fs::remove(temporary);
      throw std::runtime_error("Decomposition audit accounting failed");
    }
    fs::rename(temporary, output_path);
    const double elapsed = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    std::cout << "{\n"
              << "  \"status\": \"complete_length54_m13_native_decomposition_audit_v1\",\n"
              << "  \"representative_count\": " << expected_records << ",\n"
              << "  \"indecomposable_count\": " << indecomposable_count << ",\n"
              << "  \"decomposable_count\": " << decomposable_count << ",\n"
              << "  \"method_disagreement_count\": " << disagreement_count << ",\n"
              << "  \"elapsed_seconds\": " << elapsed << "\n"
              << "}\n";
    return 0;
  } catch (const std::exception& error) {
    if (!temporary.empty()) fs::remove(temporary);
    std::cerr << error.what() << '\n';
    return 1;
  }
}
