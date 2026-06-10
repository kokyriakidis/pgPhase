#include <iostream>
#include <string>

#include "build_catalog.hpp"
#include "collect_pipeline.hpp"
#include "graph_collect.hpp"
#include "hybrid_collect.hpp"

namespace {

void print_main_help() {
    std::cerr << "Usage: pgphase <command> [options]\n"
              << "Commands:\n"
              << "  collect-bam-variation    Collect SNP/indel evidence from a BAM/CRAM\n"
              << "  collect-graph-variation  Fast candidate collection from deconstruct VCF + GAF\n"
              << "  collect-hybrid-variation BAM calling + graph read augmentation for max phasing\n"
              << "  build-snarl-catalog      Preprocess GBZ graph into a phasing site catalog\n";
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "collect-bam-variation") {
        return collect_bam_variation(argc - 1, argv + 1);
    }
    if (argc > 1 && std::string(argv[1]) == "collect-graph-variation") {
        return collect_graph_variation(argc - 1, argv + 1);
    }
    if (argc > 1 && std::string(argv[1]) == "collect-hybrid-variation") {
        return collect_hybrid_variation(argc - 1, argv + 1);
    }
    if (argc > 1 && std::string(argv[1]) == "build-snarl-catalog") {
        return build_snarl_catalog(argc - 1, argv + 1);
    }

    print_main_help();
    return 1;
}
