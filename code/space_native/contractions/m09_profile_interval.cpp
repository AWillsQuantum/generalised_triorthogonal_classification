// Reuse the validated mathematical profile with a bounded-memory input interface.
#define main complete_m09_source_profile_main
#include "../length54/m09_source_profile_kernel.cpp"
#undef main

int main(int argc, char** argv) {
  try {
    if (argc != 5) throw std::runtime_error("Usage: INPUT OUTPUT MULTIPLICITY FIRST_INDEX");
    fs::path input_path = argv[1], output_path = argv[2];
    int multiplicity = std::stoi(argv[3]);
    std::uint64_t first = std::stoull(argv[4]);
    if (multiplicity < 0 || multiplicity > 2 || fs::exists(output_path))
      throw std::runtime_error("Invalid multiplicity or existing output");
    const auto bytes = fs::file_size(input_path);
    if (!bytes || bytes % sizeof(MaskRecord)) throw std::runtime_error("Invalid input byte count");
    const auto count = bytes / sizeof(MaskRecord);
    if (first > UINT64_MAX-count) throw std::runtime_error("Source index overflow");
    std::ifstream input(input_path, std::ios::binary);
    std::ofstream output(output_path, std::ios::binary);
    if (!input || !output) throw std::runtime_error("Could not open profile files");
    const auto quadratic = quadratic_vectors();
    for (std::uint64_t i = 0; i < count; ++i) {
      MaskRecord mask;
      SourceRecord record;
      input.read(reinterpret_cast<char*>(&mask), sizeof(mask));
      if (!input || !build_profile(mask.mask, multiplicity, first+i, 54-2*multiplicity, quadratic, record))
        throw std::runtime_error("Invalid source support or profile");
      output.write(reinterpret_cast<const char*>(&record), sizeof(record));
      if (!output) throw std::runtime_error("Profile write failed");
    }
    output.flush();
    if (!output) throw std::runtime_error("Profile flush failed");
    std::cout << "{\"source_count\":" << count << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
