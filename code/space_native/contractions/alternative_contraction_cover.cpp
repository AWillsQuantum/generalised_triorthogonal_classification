#include <algorithm>
#include <array>
#include <bit>
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
using Points = std::vector<int>;
using Histogram = std::array<int, 65>;

struct Record {
  std::array<std::uint64_t, 4> support;
  std::uint64_t orbit_size;
};
static_assert(sizeof(Record) == 40);
static_assert(std::endian::native == std::endian::little);

Points unpack(const std::array<std::uint64_t, 4>& words) {
  Points result;
  for (int i = 0; i < 4; ++i) {
    auto word = words[i];
    while (word) {
      result.push_back(64 * i + std::countr_zero(word));
      word &= word - 1;
    }
  }
  return result;
}

template <int Size>
std::array<int, Size> differences(const Points& points) {
  std::array<int, Size> result{};
  for (std::size_t i = 0; i < points.size(); ++i)
    for (std::size_t j = 0; j < i; ++j)
      result[points[i] ^ points[j]] += 2;
  return result;
}

int affine_rank(const Points& points) {
  if (points.empty()) return -1;
  std::array<int, 8> basis{};
  int rank = 0;
  for (int point : points) {
    int v = point ^ points[0];
    while (v) {
      int pivot = std::bit_width(static_cast<unsigned>(v)) - 1;
      if (basis[pivot]) v ^= basis[pivot];
      else { basis[pivot] = v; ++rank; break; }
    }
  }
  return rank;
}

bool even_moments(const Points& points, int dimension) {
  if (points.size() > 64 || points.size() % 2) return false;
  std::array<std::uint64_t, 8> rows{};
  for (std::size_t column = 0; column < points.size(); ++column)
    for (int i = 0; i < dimension; ++i)
      if (points[column] >> i & 1) rows[i] |= 1ULL << column;
  for (int i = 0; i < dimension; ++i) {
    if (std::popcount(rows[i]) % 2) return false;
    for (int j = 0; j < i; ++j) {
      if (std::popcount(rows[i] & rows[j]) % 2) return false;
      for (int k = 0; k < j; ++k)
        if (std::popcount(rows[i] & rows[j] & rows[k]) % 2) return false;
    }
  }
  return true;
}

Points contract(const Points& points, int direction) {
  if (direction < 1 || direction > 255) throw std::runtime_error("Invalid direction");
  int pivot = std::bit_width(static_cast<unsigned>(direction)) - 1;
  std::array<bool, 128> odd{};
  for (int x : points) {
    if (x >> pivot & 1) x ^= direction;
    int y = (x & ((1 << pivot) - 1)) | ((x >> (pivot + 1)) << pivot);
    odd[y] = !odd[y];
  }
  Points result;
  for (int x = 0; x < 128; ++x) if (odd[x]) result.push_back(x);
  return result;
}

Histogram histogram(const std::array<int, 128>& counts) {
  Histogram result{};
  for (int v = 1; v < 128; ++v) ++result.at(counts[v]);
  return result;
}

using LocalProfile = std::vector<std::vector<int>>;
LocalProfile local_profile(const Points& points, const std::array<int, 128>& counts) {
  LocalProfile result;
  for (int p : points) {
    std::vector<int> row;
    for (int q : points) if (p != q) row.push_back(counts[p ^ q]);
    std::sort(row.begin(), row.end());
    result.push_back(std::move(row));
  }
  std::sort(result.begin(), result.end());
  return result;
}

struct Excluded {
  int weight;
  Histogram histogram;
  LocalProfile local;
};

std::vector<Excluded> read_excluded(const fs::path& path) {
  std::ifstream file(path);
  std::string magic;
  int count = 0;
  file >> magic >> count;
  if (magic != "M8_EXCLUDED_CORES_V1" || count != 4)
    throw std::runtime_error("Invalid excluded-core header");
  std::map<int, int> weights;
  std::vector<Excluded> result;
  for (int i = 0; i < count; ++i) {
    std::uint64_t low, high;
    file >> low >> high;
    if (!file) throw std::runtime_error("Truncated excluded cores");
    Points points = unpack({low, high, 0, 0});
    if ((points.size() != 48 && points.size() != 52) || affine_rank(points) != 7 || !even_moments(points, 7))
      throw std::runtime_error("Invalid excluded core");
    auto counts = differences<128>(points);
    result.push_back({static_cast<int>(points.size()), histogram(counts), local_profile(points, counts)});
    ++weights[points.size()];
  }
  std::string extra;
  if (file >> extra) throw std::runtime_error("Unexpected excluded-core data");
  if (weights != std::map<int,int>{{48,2},{52,2}}) throw std::runtime_error("Wrong excluded-core weights");
  return result;
}

bool separates_excluded(const Points& core, const std::vector<Excluded>& excluded) {
  if (core.size() == 44) return true;
  auto counts = differences<128>(core);
  auto signature = histogram(counts);
  LocalProfile local;
  bool have_local = false;
  for (const auto& item : excluded) {
    if (item.weight != static_cast<int>(core.size()) || item.histogram != signature) continue;
    if (!have_local) { local = local_profile(core, counts); have_local = true; }
    if (local == item.local) return false;
  }
  return true;
}

bool valid_cover(const Points& points, int direction, const std::vector<Excluded>& excluded) {
  auto core = contract(points, direction);
  if (core.size() != 44 && core.size() != 48 && core.size() != 52) return false;
  if (core.size() > points.size() || affine_rank(core) != 7) return false;
  if (!even_moments(core, 7)) throw std::runtime_error("Contraction moment identity failed");
  return separates_excluded(core, excluded);
}

int find_cover(const Points& points, const std::vector<Excluded>& excluded) {
  auto counts = differences<256>(points);
  int n = static_cast<int>(points.size());
  for (int multiplicity = 0; multiplicity <= (n - 44)/2; ++multiplicity) {
    int weight = n - 2* multiplicity;
    if (weight != 44 && weight != 48 && weight != 52) continue;
    for (int direction = 2; direction < 256; ++direction)
      if (counts[direction] == 2 * multiplicity && valid_cover(points, direction, excluded)) return direction;
  }
  return 0;
}

Points read_core(const fs::path& path, int length) {
  std::ifstream file(path);
  std::string magic, identifier;
  int weight, multiplicity, count, external, affine, lift, generators;
  std::uint64_t order, fibres, raw, deficient;
  file >> magic >> identifier >> weight >> multiplicity >> count >> external
       >> affine >> lift >> generators >> order >> fibres >> raw >> deficient;
  if (!file || magic != "L" + std::to_string(length) + "_M08_NATIVE_TASK_V1"
      || count != weight || count < 0 || count > 64
      || weight + 2 * multiplicity != length || affine != 8)
    throw std::runtime_error("Invalid source input");
  Points points;
  for (int i = 0; i < count; ++i) {
    int point;
    if (!(file >> point) || point < 0 || point >= 128
        || (!points.empty() && point <= points.back()))
      throw std::runtime_error("Invalid source support");
    points.push_back(point);
  }
  if (affine_rank(points) != 7 || !even_moments(points, 7))
    throw std::runtime_error("Source support fails rank or moments");
  return points;
}

int main(int argc, char** argv) {
  try {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; i += 2) {
      if (i + 1 >= argc || !args.emplace(argv[i], argv[i+1]).second)
        throw std::runtime_error("Invalid arguments");
    }
    for (auto& [key, value] : args)
      if (key != "--input" && key != "--excluded" && key != "--length" && key != "--proof" && key != "--residual" && key != "--summary" && key != "--verify-proof" && key != "--core-input")
        throw std::runtime_error("Unknown argument");
    int length = std::stoi(args.at("--length"));
    if (length != 48 && length != 50 && length != 52 && length != 54)
      throw std::runtime_error("Unsupported target length");
    auto excluded = read_excluded(args.at("--excluded"));
    auto core = read_core(args.at("--core-input"), length);
    auto bytes = fs::file_size(args.at("--input"));
    if (!bytes || bytes % sizeof(Record)) throw std::runtime_error("Invalid candidate byte count");
    bool verify = args.contains("--verify-proof");
    std::ifstream input(args.at("--input"), std::ios::binary), proof_input;
    std::ofstream proof, residual;
    if (verify) {
      if (fs::file_size(args.at("--verify-proof")) != bytes / sizeof(Record)) throw std::runtime_error("Wrong proof byte count");
      proof_input.open(args.at("--verify-proof"), std::ios::binary);
    } else {
      for (const char* key : {"--proof", "--residual"}) {
        if (fs::exists(args.at(key))) throw std::runtime_error("Refusing to overwrite evidence");
      }
      proof.open(args.at("--proof"), std::ios::binary);
      residual.open(args.at("--residual"), std::ios::binary);
      if (!proof || !residual) throw std::runtime_error("Could not open evidence files");
    }
    if (fs::exists(args.at("--summary"))) throw std::runtime_error("Summary already exists");
    std::uint64_t total = 0, covered = 0, mass = 0, covered_mass = 0;
    std::array<std::uint64_t, 256> directions{};
    auto started = std::chrono::steady_clock::now();
    Record record;
    while (input.read(reinterpret_cast<char*>(&record), sizeof(record))) {
      Points points = unpack(record.support);
      if (points.size() != static_cast<std::size_t>(length) || !record.orbit_size || affine_rank(points) != 8 || !even_moments(points, 8))
        throw std::runtime_error("Invalid generated candidate");
      if (contract(points, 1) != core)
        throw std::runtime_error("Candidate has the wrong marked contraction");
      int direction = 0;
      if (verify) {
        direction = proof_input.get();
        if (direction < 0 || (direction && !valid_cover(points, direction, excluded)))
          throw std::runtime_error("Invalid coverage witness");
      } else {
        direction = find_cover(points, excluded);
        proof.put(static_cast<char>(direction));
        if (!direction) residual.write(reinterpret_cast<const char*>(&record), sizeof(record));
      }
      ++total; mass += record.orbit_size; ++directions[direction];
      if (direction) { ++covered; covered_mass += record.orbit_size; }
    }
    if (!input.eof() || total != bytes / sizeof(Record)) throw std::runtime_error("Candidate read failed");
    if (!verify) {
      proof.flush(); residual.flush();
      if (!proof || !residual) throw std::runtime_error("Evidence write failed");
    }
    double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::ofstream summary(args.at("--summary"));
    summary << "{\n  \"status\": \"" << (covered == total ? "complete_alternative_contraction_cover" : "unresolved_alternative_contraction_cases") << "\",\n"
            << "  \"target_length\": " << length << ",\n  \"verification_mode\": " << (verify ? "true" : "false")
            << ",\n  \"candidate_count\": " << total << ",\n  \"covered_count\": " << covered
            << ",\n  \"unresolved_count\": " << total-covered << ",\n  \"raw_marked_mass\": " << mass
            << ",\n  \"covered_raw_marked_mass\": " << covered_mass << ",\n  \"wall_seconds\": " << seconds
            << ",\n  \"direction_counts\": {";
    bool first = true;
    for (int d = 0; d < 256; ++d) if (directions[d]) {
      if (!first) summary << ',';
      first = false;
      summary << '\"' << d << "\":" << directions[d];
    }
    summary << "}\n}\n";
    if (!summary) throw std::runtime_error("Summary write failed");
    std::cout << "candidates=" << total << " covered=" << covered << " unresolved=" << total-covered << " seconds=" << seconds << '\n';
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
