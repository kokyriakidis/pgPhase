#include "collect_bam_variation.hpp"

#include "collect_pipeline.hpp"

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pgphase_collect {

enum LongOption {
    kMinAltDepthOption = 1000,
    kMinAfOption,
    kMaxAfOption,
    kReadSupportOption,
    kNoisyRegMergeDisOption,
    kMinSvLenOption,
    kChunkSizeOption,
    kOntOption,
    kStrandBiasPvalOption,
    kNoisyMaxXgapsOption,
    kMaxNoisyFracOption,
    kNoisySlideWinOption,
    kDebugSiteOption,
    kInputIsListOption
};

static std::vector<std::string> load_bam_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("failed to open BAM/CRAM list: " + path);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(in, line)) {
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;
        const size_t last = line.find_last_not_of(" \t\r\n");
        files.push_back(line.substr(first, last - first + 1));
    }
    if (files.empty()) throw std::runtime_error("BAM/CRAM list is empty: " + path);
    return files;
}

void print_collect_help() {
    std::cout << "Usage: pgphase collect-bam-variation [options] <ref.fa> <input.bam|bam.list> [region ...]\n"
              << "Options:\n"
              << "  -L, --input-is-list          Treat input path as a list of BAM/CRAM files\n"
              << "  -X, --extra-bam FILE         Extra input BAM/CRAM file; may be repeated\n"
              << "  -t, --threads INT             Region worker threads [1]\n"
              << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
              << "  -B, --min-bq INT              Minimum base quality for candidate sites [10]\n"
              << "  -D, --min-depth INT           Minimum total depth for clean candidates [5]\n"
              << "      --min-alt-depth INT       Minimum alternate depth for clean candidates [2]\n"
              << "      --min-af FLOAT            Minimum allele fraction for clean het candidates [0.20]\n"
              << "      --max-af FLOAT            Maximum allele fraction for clean het candidates [0.80]\n"
              << "  -r, --region STR              Optional region; may be repeated\n"
              << "      --region-file FILE        BED file of regions to process\n"
              << "      --autosome                Process chr1-22 / 1-22 only\n"
              << "  -j, --max-var-ratio FLOAT     Skip reads above this variant/ref-span ratio [0.05]\n"
              << "      --max-noisy-frac FLOAT    Skip reads with > this fraction in noisy regions [0.5]\n"
              << "      --include-filtered        Include QC-fail and duplicate reads\n"
              << "  -o, --output FILE             Output TSV file [output.tsv]\n"
              << "  -v, --vcf-output FILE         Optional VCF output for collected candidates\n"
              << "      --read-support FILE       Per-read ref/alt observations at candidates (for phasing)\n"
              << "      --chunk-size INT          Region chunk size in bp [500000]\n"
              << "      --noisy-merge-dis INT     Max distance (bp) to merge noisy/SV windows [500]\n"
              << "      --min-sv-len INT          min_sv_len for noisy-region cgranges merge [30]\n"
              << "      --noisy-slide-win INT     Slide window (bp) for per-read noisy regions [HiFi 100 / ONT 25]\n"
              << "      --ont                     ONT mode: Fisher exact test for alt strand bias\n"
              << "      --strand-bias-pval FLOAT  max p-value for ONT strand filter [0.01]\n"
              << "      --noisy-max-xgaps INT     max indel len (bp) for STR/homopolymer flags [5]\n"
              << "\n"
              << "Regions can also be supplied after <input.bam|bam.list>, e.g.\n"
              << "  pgphase collect-bam-variation ref.fa hifi.bam chr11:1000-2000 chr12:1-500\n";
}

} // namespace pgphase_collect

int collect_bam_variation(int argc, char* argv[]) {
    using namespace pgphase_collect;
    Options opts;
    std::vector<std::string> extra_bam_files;
    optind = 1;

    const struct option long_options[] = {
        {"threads", required_argument, nullptr, 't'},
        {"min-mapq", required_argument, nullptr, 'q'},
        {"min-bq", required_argument, nullptr, 'B'},
        {"min-depth", required_argument, nullptr, 'D'},
        {"min-alt-depth", required_argument, nullptr, kMinAltDepthOption},
        {"min-af", required_argument, nullptr, kMinAfOption},
        {"max-af", required_argument, nullptr, kMaxAfOption},
        {"region", required_argument, nullptr, 'r'},
        {"region-file", required_argument, nullptr, 'R'},
        {"autosome", no_argument, nullptr, 'a'},
        {"max-var-ratio", required_argument, nullptr, 'j'},
        {"max-noisy-frac", required_argument, nullptr, kMaxNoisyFracOption},
        {"include-filtered", no_argument, nullptr, 'f'},
        {"output", required_argument, nullptr, 'o'},
        {"vcf-output", required_argument, nullptr, 'v'},
        {"read-support", required_argument, nullptr, kReadSupportOption},
        {"chunk-size", required_argument, nullptr, kChunkSizeOption},
        {"noisy-merge-dis", required_argument, nullptr, kNoisyRegMergeDisOption},
        {"min-sv-len", required_argument, nullptr, kMinSvLenOption},
        {"noisy-slide-win", required_argument, nullptr, kNoisySlideWinOption},
        {"debug-site", required_argument, nullptr, kDebugSiteOption},
        {"extra-bam", required_argument, nullptr, 'X'},
        {"input-is-list", no_argument, nullptr, 'L'},
        {"ont", no_argument, nullptr, kOntOption},
        {"strand-bias-pval", required_argument, nullptr, kStrandBiasPvalOption},
        {"noisy-max-xgaps", required_argument, nullptr, kNoisyMaxXgapsOption},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt = 0;
    int long_index = 0;
    while ((opt = getopt_long(argc, argv, "t:q:B:D:r:R:aj:o:v:hX:L", long_options, &long_index)) != -1) {
        switch (opt) {
            case 't':
                opts.threads = std::stoi(optarg);
                break;
            case 'q':
                opts.min_mapq = std::stoi(optarg);
                break;
            case 'B':
                opts.min_bq = std::stoi(optarg);
                break;
            case 'D':
                opts.min_depth = std::stoi(optarg);
                break;
            case kMinAltDepthOption:
                opts.min_alt_depth = std::stoi(optarg);
                break;
            case kMinAfOption:
                opts.min_af = std::stod(optarg);
                break;
            case kMaxAfOption:
                opts.max_af = std::stod(optarg);
                break;
            case 'r':
                opts.regions.push_back(optarg);
                break;
            case 'R':
                opts.region_file = optarg;
                break;
            case 'a':
                opts.autosome = true;
                break;
            case 'j':
                opts.max_var_ratio_per_read = std::stod(optarg);
                break;
            case kMaxNoisyFracOption:
                opts.max_noisy_frac_per_read = std::stod(optarg);
                break;
            case 'f':
                opts.include_filtered = true;
                break;
            case 'o':
                opts.output_tsv = optarg;
                break;
            case 'v':
                opts.output_vcf = optarg;
                break;
            case kReadSupportOption:
                opts.read_support_tsv = optarg;
                break;
            case kChunkSizeOption:
                opts.chunk_size = std::stoll(optarg);
                break;
            case kNoisyRegMergeDisOption:
                opts.noisy_reg_merge_dis = std::stoi(optarg);
                break;
            case kMinSvLenOption:
                opts.min_sv_len = std::stoi(optarg);
                break;
            case kNoisySlideWinOption:
                opts.noisy_reg_slide_win = std::stoi(optarg);
                break;
            case kDebugSiteOption:
                opts.debug_site = optarg;
                break;
            case 'X':
                extra_bam_files.push_back(optarg);
                break;
            case 'L':
                opts.input_is_list = true;
                break;
            case kOntOption:
                opts.is_ont = true;
                break;
            case kStrandBiasPvalOption:
                opts.strand_bias_pval = std::stod(optarg);
                break;
            case kNoisyMaxXgapsOption:
                opts.noisy_reg_max_xgaps = std::stoi(optarg);
                break;
            case 'h':
                print_collect_help();
                return 0;
            default:
                print_collect_help();
                return 1;
        }
    }

    if (opts.threads < 1 || opts.min_mapq < 0 || opts.min_bq < 0 || opts.chunk_size < 1 ||
        opts.min_depth < 0 || opts.min_alt_depth < 0 || opts.noisy_reg_merge_dis < 0 || opts.min_sv_len < 0 ||
        opts.min_af < 0.0 || opts.max_af < opts.min_af || opts.strand_bias_pval < 0.0 ||
        opts.strand_bias_pval > 1.0 || opts.max_var_ratio_per_read < 0.0 || opts.max_noisy_frac_per_read < 0.0 ||
        opts.noisy_reg_max_xgaps < 0 || opts.noisy_reg_slide_win < -1) {
        std::cerr << "Error: numeric thresholds are invalid\n";
        return 1;
    }

    if (optind + 2 > argc) {
        print_collect_help();
        return 1;
    }
    opts.ref_fasta = argv[optind];
    const std::string input_path = argv[optind + 1];
    if (opts.input_is_list) {
        opts.bam_files = load_bam_list(input_path);
    } else {
        opts.bam_files.push_back(input_path);
    }
    opts.bam_files.insert(opts.bam_files.end(), extra_bam_files.begin(), extra_bam_files.end());
    opts.bam_file = opts.bam_files.front();
    for (int arg_i = optind + 2; arg_i < argc; ++arg_i) {
        opts.regions.push_back(argv[arg_i]);
    }

    try {
        run_collect_bam_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
