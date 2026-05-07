#ifndef PGPHASE_PHASE_ENGINE_HPP
#define PGPHASE_PHASE_ENGINE_HPP

#include "collect_types.hpp"

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

namespace pgphase_collect {

constexpr int kPhaseAlleleMissing = -1;
constexpr int kPhaseAlleleAmbiguous = -2;

struct PhaseObservation {
    int site_index = -1;
    int allele = kPhaseAlleleMissing;
};

struct PhaseSite {
    std::string site_id;
    std::string chrom;
    hts_pos_t order_pos = 0;
    int n_alleles = 0;
    std::vector<int> allele_counts;
    hts_pos_t phase_set = -1;
    std::array<int, 3> hap_to_cons_allele{-1, -1, -1};
    int score = 1;
    bool eligible = true;

    std::string parent_site_id;
    std::string root_site_id;
    int level = -1;
    std::vector<int> conditional_parent_alleles;
};

struct ReadPhaseProfile {
    std::string read_name;
    std::vector<PhaseObservation> observations;
    int hap = 0;
    hts_pos_t phase_set = -1;
};

struct PhaseMatrix {
    std::vector<PhaseSite> sites;
    std::vector<ReadPhaseProfile> reads;
};

struct PhaseChunk {
    std::string chrom;
    hts_pos_t beg = 0; // 0-based half-open.
    hts_pos_t end = 0;
    PhaseMatrix matrix;
};

void phase_matrix(PhaseMatrix& matrix, int max_iters = 10);

void stitch_phase_chunks(std::vector<PhaseChunk>& chunks);

void write_phase_site_counts_tsv(std::ostream& out, const PhaseMatrix& matrix);

void write_phase_read_profiles_tsv(std::ostream& out, const PhaseMatrix& matrix);

void write_phase_sites_tsv(std::ostream& out, const PhaseMatrix& matrix);

void write_phase_reads_tsv(std::ostream& out, const PhaseMatrix& matrix);

} // namespace pgphase_collect

#endif
