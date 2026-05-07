#include "build_catalog.hpp"

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

namespace {

std::string bc_shell_quote(const std::string& val) {
    std::string out = "'";
    for (char c : val) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

void print_help(const char* prog) {
    std::cerr
        << "Usage: " << prog << " build-snarl-catalog [options] <graph.gbz>\n"
        << "\n"
        << "Run 'vg deconstruct -a' on a GBZ pangenome graph and write a bgzipped,\n"
        << "tabix-indexed VCF that serves as the phasing site catalog.\n"
        << "Use the output directly with 'phase-graph --graph-sites-vcf'.\n"
        << "\n"
        << "  <graph.gbz> may be the full pangenome GBZ or a per-chr chunk extracted\n"
        << "  with 'vg chunk --gbz --contig CHR -x full.gbz -o /tmp/chunk'.\n"
        << "  Using a per-chr chunk keeps 'vg deconstruct' RAM under ~1 GB per chromosome.\n"
        << "\n"
        << "Options:\n"
        << "  --ref-sample STR    Reference sample name, e.g. CHM13 [required]\n"
        << "  --contig STR        Process one contig, e.g. chr20 (recommended: per-chr)\n"
        << "  -t, --threads N     Threads for vg deconstruct [4]\n"
        << "  --vg-bin PATH       Path to vg binary [vg]\n"
        << "  --snarls FILE       Pre-built snarls file from 'vg snarls' (skips recomputation)\n"
        << "  -o, --output FILE   Output bgzipped VCF file (tabix index created automatically) [required]\n"
        << "  -h, --help          Show this help\n"
        << "\n"
        << "Example — per-chr (low RAM, recommended):\n"
        << "  # Precompute snarls once\n"
        << "  vg snarls -t 16 full.gbz > full.snarls.pb\n"
        << "\n"
        << "  # Extract per-chr chunks and build catalog VCF in parallel\n"
        << "  for chr in chr{1..22} chrX chrY; do\n"
        << "    vg chunk --gbz --contig $chr -x full.gbz -o /tmp/chunk_${chr}\n"
        << "    pgphase build-snarl-catalog --ref-sample CHM13 --contig $chr \\\n"
        << "        --snarls full.snarls.pb --threads 8 \\\n"
        << "        -o ${chr}.sites.vcf.gz \\\n"
        << "        /tmp/chunk_${chr}_graph_0_${chr}.gbz &\n"
        << "  done\n"
        << "  wait\n"
        << "\n"
        << "  # Run phasing\n"
        << "  pgphase phase-graph --graph-sites-vcf chr20.sites.vcf.gz ...\n";
}

} // namespace

int build_snarl_catalog(int argc, char* argv[]) {
    std::string ref_sample;
    std::string contig;
    int threads = 4;
    std::string vg_bin = "vg";
    std::string snarls_file;
    std::string output_file;

    static const option long_opts[] = {
        {"ref-sample", required_argument, nullptr, 'R'},
        {"contig",     required_argument, nullptr, 'N'},
        {"threads",    required_argument, nullptr, 't'},
        {"vg-bin",     required_argument, nullptr, 'V'},
        {"snarls",     required_argument, nullptr, 'S'},
        {"output",     required_argument, nullptr, 'o'},
        {"help",       no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int ch;
    while ((ch = getopt_long(argc, argv, "t:o:h", long_opts, nullptr)) != -1) {
        switch (ch) {
        case 'R': ref_sample  = optarg; break;
        case 'N': contig      = optarg; break;
        case 't': threads     = std::stoi(optarg); break;
        case 'V': vg_bin      = optarg; break;
        case 'S': snarls_file = optarg; break;
        case 'o': output_file = optarg; break;
        case 'h': print_help(argv[0]); return 0;
        default:  print_help(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        std::cerr << "error: missing <graph.gbz>\n";
        print_help(argv[0]);
        return 1;
    }
    if (ref_sample.empty()) {
        std::cerr << "error: --ref-sample is required (e.g. --ref-sample CHM13)\n";
        return 1;
    }
    if (output_file.empty()) {
        std::cerr << "error: -o/--output is required\n";
        print_help(argv[0]);
        return 1;
    }

    const std::string gbz_file = argv[optind];

    // Build the vg deconstruct command.
    std::ostringstream deconstruct_cmd;
    deconstruct_cmd << bc_shell_quote(vg_bin) << " deconstruct";
    if (!contig.empty())
        deconstruct_cmd << " -p " << bc_shell_quote(ref_sample + "#0#" + contig);
    else
        deconstruct_cmd << " -P " << bc_shell_quote(ref_sample + "#0#");
    deconstruct_cmd << " -a -t " << threads;
    if (!snarls_file.empty())
        deconstruct_cmd << " -r " << bc_shell_quote(snarls_file);
    deconstruct_cmd << " " << bc_shell_quote(gbz_file);

    const std::string full_cmd =
        deconstruct_cmd.str() + " | bgzip -c > " + bc_shell_quote(output_file);
    std::cerr << "[build-snarl-catalog] " << full_cmd << "\n";
    int ret = std::system(full_cmd.c_str());

    if (ret == 0) {
        const std::string tabix_cmd = "tabix -f -p vcf " + bc_shell_quote(output_file);
        std::cerr << "[build-snarl-catalog] " << tabix_cmd << "\n";
        ret = std::system(tabix_cmd.c_str());
    }

    if (ret != 0)
        std::cerr << "[build-snarl-catalog] warning: command exited with non-zero status\n";
    return ret == 0 ? 0 : 1;
}
