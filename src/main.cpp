#include <iostream>
#include <string>

#include "collect_pipeline.hpp"

namespace {

void print_main_help() {
    std::cerr << "Usage: pgphase <command> [options]\n"
              << "Commands:\n"
              << "  collect-bam-variation    Collect SNP/indel evidence from a BAM/CRAM\n";
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "collect-bam-variation") {
        return collect_bam_variation(argc - 1, argv + 1);
    }

    print_main_help();
    return 1;
}
