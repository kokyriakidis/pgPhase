/**
 * @file hybrid_collect.cpp
 * @brief CLI entry point for `collect-hybrid-variation`.
 *
 * Parses options and delegates to run_collect_hybrid_variation.
 * Accepts all BAM pipeline options plus --graph-sites and --gaf.
 */

#include "hybrid_collect.hpp"
#include "collect_pipeline.hpp"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <string>

namespace pgphase_collect {

static void print_hybrid_help() {
    std::cout
        << "Usage: pgphase collect-hybrid-variation [options]\n"
        << "\n"
        << "Hybrid BAM+graph phasing: BAM variant calling augmented with graph\n"
        << "snarl site observations for maximum read phasing coverage.\n"
        << "\n"
        << "Required:\n"
        << "      --ref FILE                Reference FASTA (indexed)\n"
        << "      --bam FILE                Input BAM/CRAM file\n"
        << "      --graph-sites FILE        Graph sites VCF (tabix-indexed, from vg deconstruct)\n"
        << "      --gaf FILE                Coordinate-sorted GAF (bgzipped + tabix-indexed)\n"
        << "\n"
        << "Options:\n"
        << "  -t, --threads INT             Region worker threads [1]\n"
        << "  -q, --min-mapq INT            Minimum read mapping quality [30]\n"
        << "  -B, --min-bq INT              Minimum base quality [10]\n"
        << "  -D, --min-depth INT           Minimum total depth [5]\n"
        << "      --min-alt-depth INT       Minimum alternate depth [2]\n"
        << "      --min-af FLOAT            Minimum allele fraction [0.20]\n"
        << "      --max-af FLOAT            Maximum allele fraction [0.80]\n"
        << "  -r, --region STR              Region to process; may be repeated\n"
        << "      --region-file FILE        BED file of regions\n"
        << "      --autosome                Process chr1-22 / 1-22 only\n"
        << "      --chunk-size INT          Region chunk size in bp [500000]\n"
        << "      --hifi                    HiFi mode [default]\n"
        << "      --ont                     ONT mode\n"
        << "  -o, --output FILE             Output TSV [output.tsv]\n"
        << "  -v, --vcf-output FILE         Candidate VCF output\n"
        << "      --phased-vcf-out FILE     Phased VCF (GT:DP:AD:VAF:GQ:PS)\n"
        << "  -b, --out-bam FILE            Phased BAM output\n"
        << "      --pgbam-file FILE         Optional .pgbam sidecar for stitching\n"
        << "      --stitch-min-margin INT   Min flip-vote margin to merge blocks [10]\n"
        << "      --stitch-rule INT         Stitch rule: 0=net-margin 1=both-strands 2=literal 3=both+margin [1]\n"
        << "      --graph-indel-af-margin F Max |AF-0.5| for graph het-indel anchor [0.11]\n"
        << "      --graph-indel-min-alt INT Min alt support for graph het-indel anchor [0]\n"
        << "      --keep-noisy-kmeans       Restore step-4 noisy-candidate k-means re-orientation (off by default)\n"
        << "      --no-hybrid-trim          Disable minimal-VCF trimming of graph-only alleles before noise filter (on by default)\n"
        << "      --gap-fill                Additively phase reads the clean core left unphased into a disjoint PS namespace (off by default)\n"
        << "  -V, --verbose INT             Verbosity level [0]\n"
        << "\n"
        << "Examples:\n"
        << "  pgphase collect-hybrid-variation \\\n"
        << "      --ref grch38.fa \\\n"
        << "      --bam surjected.bam \\\n"
        << "      --graph-sites sites.vcf.gz \\\n"
        << "      --gaf reads.gaf.gz \\\n"
        << "      --phased-vcf-out phased.vcf \\\n"
        << "      -b phased.bam \\\n"
        << "      -t 8\n";
}

} // namespace pgphase_collect

int pgphase_collect::collect_hybrid_variation(int argc, char* argv[]) {
    Options opts;
    // Hybrid uses a non-zero stitch margin by default (see
    // kHybridDefaultStitchMinMargin); --stitch-min-margin overrides it.
    opts.stitch_min_margin = kHybridDefaultStitchMinMargin;
    // Hybrid defaults to the both-strands-bridged stitch rule (see
    // kHybridDefaultStitchRule); --stitch-rule overrides it.
    opts.stitch_rule = kHybridDefaultStitchRule;
    // Hybrid skips the step-4 noisy-candidate k-means re-orientation by default:
    // it phased ~8k extra reads at ~65% error and poisoned the BAM-shared core.
    // --keep-noisy-kmeans restores the old behaviour. See CHECKPOINT.md.
    opts.skip_noisy_kmeans = true;
    // Hybrid trims graph-only catalog alleles to minimal VCF form before the
    // indel noise filter by default, matching apply_graph_noise_filter. Catalog
    // alleles carry the full repeat run on both flanks, so the untrimmed filter
    // over-demoted ~342 real het indels to REP_HET_INDEL (excluded from k-means);
    // trimming recovers them as anchors and lowers hamming/switch error.
    // --no-hybrid-trim restores the untrimmed behaviour. See CHECKPOINT.md.
    opts.exp_hybrid_trim = true;

    enum HybridLongOption {
        kRefOption = 1000,
        kBamOption,
        kGraphSitesOption,
        kGafOption,
        kMinAltDepthOption,
        kMinAfOption,
        kMaxAfOption,
        kChunkSizeOption,
        kHifiOption,
        kOntOption,
        kPhasedVcfOutputOption,
        kPgbamFileOption,
        kAutosomeOption,
        kRegionFileOption,
        kNoisyMaxXgapsOption,
        kRefineAlnOption,
        kStitchMinMarginOption,
        kStitchRuleOption,
        kGraphIndelAfMarginOption,
        kGraphIndelMinAltOption,
        kKeepNoisyKmeansOption,
        kNoHybridTrimOption,
        kGapFillOption,
    };

    static struct option long_options[] = {
        {"ref",             required_argument, nullptr, kRefOption},
        {"bam",             required_argument, nullptr, kBamOption},
        {"graph-sites",     required_argument, nullptr, kGraphSitesOption},
        {"gaf",             required_argument, nullptr, kGafOption},
        {"threads",         required_argument, nullptr, 't'},
        {"min-mapq",        required_argument, nullptr, 'q'},
        {"min-bq",          required_argument, nullptr, 'B'},
        {"min-depth",       required_argument, nullptr, 'D'},
        {"min-alt-depth",   required_argument, nullptr, kMinAltDepthOption},
        {"min-af",          required_argument, nullptr, kMinAfOption},
        {"max-af",          required_argument, nullptr, kMaxAfOption},
        {"region",          required_argument, nullptr, 'r'},
        {"region-file",     required_argument, nullptr, kRegionFileOption},
        {"autosome",        no_argument,       nullptr, kAutosomeOption},
        {"chunk-size",      required_argument, nullptr, kChunkSizeOption},
        {"hifi",            no_argument,       nullptr, kHifiOption},
        {"ont",             no_argument,       nullptr, kOntOption},
        {"output",          required_argument, nullptr, 'o'},
        {"vcf-output",      required_argument, nullptr, 'v'},
        {"phased-vcf-out",  required_argument, nullptr, kPhasedVcfOutputOption},
        {"out-bam",         required_argument, nullptr, 'b'},
        {"pgbam-file",      required_argument, nullptr, kPgbamFileOption},
        {"noisy-max-xgaps", required_argument, nullptr, kNoisyMaxXgapsOption},
        {"stitch-min-margin", required_argument, nullptr, kStitchMinMarginOption},
        {"stitch-rule", required_argument, nullptr, kStitchRuleOption},
        {"graph-indel-af-margin", required_argument, nullptr, kGraphIndelAfMarginOption},
        {"graph-indel-min-alt", required_argument, nullptr, kGraphIndelMinAltOption},
        {"keep-noisy-kmeans", no_argument,     nullptr, kKeepNoisyKmeansOption},
        {"no-hybrid-trim", no_argument,        nullptr, kNoHybridTrimOption},
        {"gap-fill",        no_argument,       nullptr, kGapFillOption},
        {"refine-aln",      no_argument,       nullptr, kRefineAlnOption},
        {"verbose",         required_argument, nullptr, 'V'},
        {"help",            no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "t:q:B:D:r:o:v:b:V:h",
                            long_options, nullptr)) != -1) {
        switch (c) {
            case kRefOption:          opts.ref_fasta = optarg; break;
            case kBamOption:          opts.bam_files.push_back(optarg); break;
            case kGraphSitesOption:   opts.graph_sites_vcf = optarg; break;
            case kGafOption:          opts.gaf_file = optarg; break;
            case 't':                 opts.threads = std::atoi(optarg); break;
            case 'q':                 opts.min_mapq = std::atoi(optarg); break;
            case 'B':                 opts.min_bq = std::atoi(optarg); break;
            case 'D':                 opts.min_depth = std::atoi(optarg); break;
            case kMinAltDepthOption:  opts.min_alt_depth = std::atoi(optarg); break;
            case kMinAfOption:        opts.min_af = std::atof(optarg); break;
            case kMaxAfOption:        opts.max_af = std::atof(optarg); break;
            case 'r':                 opts.regions.push_back(optarg); break;
            case kRegionFileOption:   opts.region_file = optarg; break;
            case kAutosomeOption:     opts.autosome = true; break;
            case kChunkSizeOption:    opts.chunk_size = std::atoi(optarg); break;
            case kHifiOption:
                opts.read_technology = ReadTechnology::Hifi;
                break;
            case kOntOption:
                opts.read_technology = ReadTechnology::Ont;
                break;
            case 'o':                 opts.output_tsv = optarg; break;
            case 'v':                 opts.output_vcf = optarg; break;
            case kPhasedVcfOutputOption: opts.output_phased_vcf = optarg; break;
            case 'b':
                opts.output_aln = optarg;
                opts.output_aln_format = OutputAlignmentFormat::Bam;
                break;
            case kPgbamFileOption:    opts.pgbam_file = optarg; break;
            case kNoisyMaxXgapsOption: opts.noisy_reg_max_xgaps = std::atoi(optarg); break;
            case kStitchMinMarginOption: opts.stitch_min_margin = std::atoi(optarg); break;
            case kStitchRuleOption: opts.stitch_rule = std::atoi(optarg); break;
            case kGraphIndelAfMarginOption: opts.graph_indel_af_margin = std::atof(optarg); break;
            case kGraphIndelMinAltOption: opts.graph_indel_min_alt = std::atoi(optarg); break;
            case kKeepNoisyKmeansOption: opts.skip_noisy_kmeans = false; break;
            case kNoHybridTrimOption: opts.exp_hybrid_trim = false; break;
            case kGapFillOption:      opts.gap_fill = true; break;
            case kRefineAlnOption:    opts.refine_aln = true; break;
            case 'V':                 opts.verbose = std::atoi(optarg); break;
            case 'h':
                print_hybrid_help();
                return 0;
            default:
                print_hybrid_help();
                return 1;
        }
    }

    if (opts.ref_fasta.empty() || opts.bam_files.empty() ||
        opts.graph_sites_vcf.empty() || opts.gaf_file.empty()) {
        std::cerr << "Error: --ref, --bam, --graph-sites, and --gaf are required\n\n";
        print_hybrid_help();
        return 1;
    }

    try {
        run_collect_hybrid_variation(opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
