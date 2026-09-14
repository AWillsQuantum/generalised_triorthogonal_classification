// Exact affine equivalence via the point/hyperplane incidence structure.
#include <bliss/graph.hh>
#include <bliss/stats.hh>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef BLISS_USE_GMP
#error Exact stabiliser orders require BLISS_USE_GMP
#endif

using Points = std::vector<unsigned>;

unsigned rank(const Points& vectors, unsigned m) {
  Points rows(m);
  unsigned result = 0;
  for (auto value : vectors) {
    while (value) {
      const auto pivot = std::bit_width(value) - 1;
      if (rows[pivot]) value ^= rows[pivot];
      else { rows[pivot] = value; ++result; break; }
    }
  }
  return result;
}

void array(const Points& values) {
  std::cout << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << values[i];
  }
  std::cout << ']';
}

Points support_bits(std::string text, unsigned size) {
  if (text.starts_with("0x")) text.erase(0, 2);
  if (text.empty()) throw std::runtime_error("Empty truth table");
  Points result(size);
  unsigned position = 0;
  for (auto p = text.rbegin(); p != text.rend(); ++p) {
    unsigned digit;
    if (*p >= '0' && *p <= '9') digit = *p - '0';
    else if (*p >= 'a' && *p <= 'f') digit = *p - 'a' + 10;
    else if (*p >= 'A' && *p <= 'F') digit = *p - 'A' + 10;
    else throw std::runtime_error("Invalid hexadecimal truth table");
    for (unsigned bit = 0; bit < 4; ++bit, ++position) {
      if (position < size) result[position] = (digit >> bit) & 1;
      else if ((digit >> bit) & 1) throw std::runtime_error("Truth table exceeds dimension");
    }
  }
  return result;
}

void canonicalise(unsigned m, const std::string& mask) {
  if (m < 2 || m > 10) throw std::runtime_error("Graph backend requires 2 <= m <= 10");
  const unsigned size = 1U << m;
  const auto bits = support_bits(mask, size);
  bliss::Graph graph;
  for (unsigned x = 0; x < size; ++x) graph.add_vertex(bits[x]);
  for (unsigned a = 1; a < size; ++a) {
    for (unsigned b = 0; b < 2; ++b) {
      const auto h = graph.add_vertex(2);
      for (unsigned x = 0; x < size; ++x)
        if ((std::popcount(a & x) & 1U) == b) graph.add_edge(x, h);
    }
  }
  graph.set_splitting_heuristic(bliss::Graph::shs_fsm);
  std::vector<Points> generators;
  auto report = [&](unsigned n, const unsigned* permutation) {
    if (n != 3 * size - 2) throw std::runtime_error("Incorrect graph order");
    const auto origin = permutation[0];
    Points basis;
    for (unsigned i = 0; i < m; ++i) basis.push_back(permutation[1U << i] ^ origin);
    if (origin >= size || rank(basis, m) != m)
      throw std::runtime_error("Non-affine graph generator");
    for (unsigned x = 0; x < size; ++x) {
      unsigned y = origin;
      for (unsigned i = 0; i < m; ++i) if ((x >> i) & 1) y ^= basis[i];
      if (y >= size || permutation[x] != y || bits[x] != bits[y])
        throw std::runtime_error("Generator fails affine support verification");
    }
    Points frame{origin};
    frame.insert(frame.end(), basis.begin(), basis.end());
    generators.push_back(frame);
  };
  bliss::Stats stats;
  const auto* labels = graph.canonical_form(stats, report);
  Points point_order(size);
  std::iota(point_order.begin(), point_order.end(), 0);
  std::sort(point_order.begin(), point_order.end(),
            [&](unsigned a, unsigned b) { return labels[a] < labels[b]; });
  const auto origin = point_order[0];
  Points basis;
  for (auto point : point_order) {
    auto extended = basis;
    extended.push_back(point ^ origin);
    if (rank(extended, m) > basis.size()) basis = extended;
    if (basis.size() == m) break;
  }
  if (basis.size() != m) throw std::runtime_error("Incomplete canonical affine frame");
  std::string canonical((size + 3) / 4, '0');
  constexpr char digits[] = "0123456789abcdef";
  Points nibbles((size + 3) / 4);
  for (unsigned x = 0; x < size; ++x) {
    unsigned y = origin;
    for (unsigned i = 0; i < m; ++i) if ((x >> i) & 1) y ^= basis[i];
    nibbles[x / 4] |= bits[y] << (x % 4);
  }
  for (unsigned i = 0; i < nibbles.size(); ++i)
    canonical[canonical.size() - 1 - i] = digits[nibbles[i]];
  mpz_t group;
  mpz_init(group);
  stats.get_group_size().get(group);
  std::vector<char> group_text(mpz_sizeinbase(group, 10) + 2);
  mpz_get_str(group_text.data(), 10, group);
  mpz_clear(group);
  Points frame{origin};
  frame.insert(frame.end(), basis.begin(), basis.end());
  std::cout << "{\"m\":" << m << ",\"canonical_mask\":\"" << canonical
            << "\",\"stabilizer_order\":" << group_text.data() << ",\"frame\":";
  array(frame);
  std::cout << ",\"generators\":[";
  for (std::size_t i = 0; i < generators.size(); ++i) {
    if (i) std::cout << ',';
    array(generators[i]);
  }
  std::cout << "]}\n";
}

int main() {
  try {
    unsigned m;
    std::string mask;
    while (std::cin >> m) {
      if (!(std::cin >> mask)) throw std::runtime_error("Truncated input record");
      canonicalise(m, mask);
    }
    if (!std::cin.eof()) throw std::runtime_error("Invalid input dimension");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
