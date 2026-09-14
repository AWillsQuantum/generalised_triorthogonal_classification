#define main quotient_kernel_main
#include "../../code/space_native/length54/m13_exact_quotient_kernel.cpp"
#undef main

int main() {
  try {
    int first;
    while (std::cin >> first) {
      Mask mask;
      for (int i = 0; i < kSupportSize; ++i) {
        int point = first;
        if (i && !(std::cin >> point)) throw std::runtime_error("Truncated support");
        if (point < 0 || point >= kAmbientSize || point_is_set(mask, point)) {
          throw std::runtime_error("Invalid support point");
        }
        mask.words[point / 64] |= 1ULL << (point % 64);
      }
      const auto data = prepare_support(mask);
      std::cout << '[';
      for (int i = 0; i < kSupportSize; ++i) {
        if (i) std::cout << ',';
        std::cout << static_cast<int>(data.colours[i]);
      }
      std::cout << "]\n";
    }
    if (!std::cin.eof()) throw std::runtime_error("Invalid support input");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
