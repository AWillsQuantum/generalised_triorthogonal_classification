// Exact affine equivalence from a spanning set of low-weight codewords.
#include <bliss/graph.hh>
#include <bliss/stats.hh>
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Word = std::uint64_t;
using Points = std::vector<unsigned>;

bool insert(Word value, std::array<Word, 64>& basis) {
  while (value) {
    const auto pivot = std::bit_width(value)-1;
    if (basis[pivot]) value ^= basis[pivot];
    else { basis[pivot] = value; return true; }
  }
  return false;
}

void array(const Points& values) {
  std::cout << '[';
  for (unsigned i = 0; i < values.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << values[i];
  }
  std::cout << ']';
}

void canonicalise(unsigned m, const Points& points) {
  const auto n = points.size();
  if (m < 1 || m > 17 || n < m+1 || n > 63 || !std::is_sorted(points.begin(), points.end()) ||
      std::adjacent_find(points.begin(), points.end()) != points.end() || points.back() >= (1U << m)) {
    throw std::runtime_error("Invalid full-rank support domain");
  }
  std::vector<Word> rows(m+1);
  rows[0] = (Word{1} << n)-1;
  for (unsigned j = 0; j < n; ++j)
    for (unsigned i = 0; i < m; ++i) rows[i+1] |= Word((points[j] >> i) & 1) << j;
  std::array<Word, 64> row_basis{};
  for (auto row : rows)
    if (!insert(row, row_basis)) throw std::runtime_error("The support lacks full affine rank");
  std::vector<Word> words;
  words.reserve((1U << (m+1))-1);
  Word word = 0;
  for (unsigned i = 1; i < (1U << (m+1)); ++i) {
    word ^= rows[std::countr_zero(i)];
    words.push_back(word);
  }
  std::sort(words.begin(), words.end(), [](Word a, Word b) {
    if (std::popcount(a) != std::popcount(b)) return std::popcount(a) < std::popcount(b);
    return a < b;
  });
  std::array<Word, 64> span{};
  unsigned rank = 0, maximum_weight = 0;
  for (auto value : words) {
    if (insert(value, span)) ++rank;
    if (rank == m+1) { maximum_weight = std::popcount(value); break; }
  }
  if (rank != m+1) throw std::runtime_error("The codewords do not span the input");
  const auto end = std::find_if(words.begin(), words.end(), [&](Word x) {
    return std::popcount(x) > maximum_weight;
  });
  words.erase(end, words.end());
  bliss::Graph graph;
  for (unsigned j = 0; j < n; ++j) graph.add_vertex(0);
  for (auto value : words) {
    const auto vertex = graph.add_vertex(1+std::popcount(value));
    while (value) {
      const auto j = std::countr_zero(value);
      graph.add_edge(j, vertex);
      value &= value-1;
    }
  }
  graph.set_splitting_heuristic(bliss::Graph::shs_fsm);
  bliss::Stats stats;
  const auto* labels = graph.canonical_form(stats, [](unsigned, const unsigned*) {});
  Points order(n);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](unsigned a, unsigned b) { return labels[a] < labels[b]; });
  const auto origin = points[order[0]];
  std::array<Word, 64> affine_span{};
  Points basis;
  for (auto i : order) {
    const auto value = points[i] ^ origin;
    if (insert(value, affine_span)) basis.push_back(value);
    if (basis.size() == m) break;
  }
  if (basis.size() != m) throw std::runtime_error("The canonical affine frame is incomplete");
  Points images(1, 0);
  for (auto value : basis) {
    const auto size = images.size();
    for (unsigned i = 0; i < size; ++i) images.push_back(images[i] ^ value);
  }
  Points inverse(1U << m);
  for (unsigned i = 0; i < images.size(); ++i) inverse[images[i]] = i;
  Points canonical;
  for (auto point : points) canonical.push_back(inverse[point ^ origin]);
  std::sort(canonical.begin(), canonical.end());
  Points frame{origin};
  frame.insert(frame.end(), basis.begin(), basis.end());
  std::cout << "{\"m\":" << m << ",\"n\":" << n << ",\"maximum_codeword_weight\":" << maximum_weight
            << ",\"selected_codewords\":" << words.size() << ",\"canonical_points\":";
  array(canonical);
  std::cout << ",\"frame\":";
  array(frame);
  std::cout << "}\n";
}

int main() {
  try {
    unsigned m, n;
    while (std::cin >> m) {
      if (!(std::cin >> n) || n > 63) throw std::runtime_error("Invalid support length");
      Points points(n);
      for (auto& point : points)
        if (!(std::cin >> point)) throw std::runtime_error("Truncated support");
      canonicalise(m, points);
    }
    if (!std::cin.eof()) throw std::runtime_error("Invalid support dimension");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
