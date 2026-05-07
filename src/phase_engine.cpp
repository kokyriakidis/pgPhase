#include "phase_engine.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <ostream>
#include <unordered_map>
#include <vector>

namespace pgphase_collect {

namespace {

bool site_is_informative(const PhaseSite& site) {
    if (!site.eligible || site.n_alleles < 2) return false;
    int nonzero = 0;
    for (int count : site.allele_counts) {
        if (count > 0) ++nonzero;
    }
    return nonzero >= 2;
}

void init_site_consensus(PhaseMatrix& matrix) {
    for (PhaseSite& site : matrix.sites) {
        site.phase_set = site_is_informative(site) ? site.order_pos : -1;
        site.hap_to_cons_allele = {-1, -1, -1};
        if (!site_is_informative(site)) continue;

        int best1 = -1;
        int best2 = -1;
        for (int allele = 0; allele < site.n_alleles; ++allele) {
            if (best1 == -1 ||
                site.allele_counts[static_cast<size_t>(allele)] >
                    site.allele_counts[static_cast<size_t>(best1)]) {
                best2 = best1;
                best1 = allele;
            } else if (best2 == -1 ||
                       site.allele_counts[static_cast<size_t>(allele)] >
                           site.allele_counts[static_cast<size_t>(best2)]) {
                best2 = allele;
            }
        }
        site.hap_to_cons_allele[1] = best1;
        site.hap_to_cons_allele[2] = best2;
    }
}

int score_read_hap(const PhaseMatrix& matrix,
                   const ReadPhaseProfile& read,
                   int hap) {
    int score = 0;
    const int other = 3 - hap;
    for (const PhaseObservation& obs : read.observations) {
        if (obs.site_index < 0 ||
            static_cast<size_t>(obs.site_index) >= matrix.sites.size()) {
            continue;
        }
        const PhaseSite& site = matrix.sites[static_cast<size_t>(obs.site_index)];
        const int hap_allele = site.hap_to_cons_allele[static_cast<size_t>(hap)];
        const int other_allele = site.hap_to_cons_allele[static_cast<size_t>(other)];
        if (hap_allele < 0 || other_allele < 0 || hap_allele == other_allele) continue;
        if (obs.allele == hap_allele) score += site.score;
        else if (obs.allele == other_allele) score -= site.score;
    }
    return score;
}

int assign_read_hap(const PhaseMatrix& matrix, const ReadPhaseProfile& read) {
    const int s1 = score_read_hap(matrix, read, 1);
    const int s2 = score_read_hap(matrix, read, 2);
    if (s1 == 0 && s2 == 0) return 0;
    if (s1 == s2) return 0;
    return (s1 > s2) ? 1 : 2;
}

hts_pos_t first_read_phase_set(const PhaseMatrix& matrix,
                               const ReadPhaseProfile& read) {
    for (const PhaseObservation& obs : read.observations) {
        if (obs.site_index < 0 ||
            static_cast<size_t>(obs.site_index) >= matrix.sites.size()) {
            continue;
        }
        const PhaseSite& site = matrix.sites[static_cast<size_t>(obs.site_index)];
        if (site.phase_set > 0 &&
            site.hap_to_cons_allele[1] >= 0 &&
            site.hap_to_cons_allele[2] >= 0 &&
            site.hap_to_cons_allele[1] != site.hap_to_cons_allele[2]) {
            return site.phase_set;
        }
    }
    return -1;
}

bool site_is_phased(const PhaseSite& site) {
    return site.phase_set > 0 &&
           site.hap_to_cons_allele[1] >= 0 &&
           site.hap_to_cons_allele[2] >= 0 &&
           site.hap_to_cons_allele[1] != site.hap_to_cons_allele[2];
}

int find_parent(std::vector<int>& parent, int x) {
    while (parent[static_cast<size_t>(x)] != x) {
        parent[static_cast<size_t>(x)] =
            parent[static_cast<size_t>(parent[static_cast<size_t>(x)])];
        x = parent[static_cast<size_t>(x)];
    }
    return x;
}

void union_parent(std::vector<int>& parent, int a, int b) {
    const int ra = find_parent(parent, a);
    const int rb = find_parent(parent, b);
    if (ra != rb) parent[static_cast<size_t>(rb)] = ra;
}

void assign_phase_sets_from_read_components(PhaseMatrix& matrix) {
    std::vector<int> parent(matrix.sites.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);

    for (const ReadPhaseProfile& read : matrix.reads) {
        if (read.hap == 0) continue;
        int first_site = -1;
        for (const PhaseObservation& obs : read.observations) {
            if (obs.site_index < 0 ||
                static_cast<size_t>(obs.site_index) >= matrix.sites.size()) {
                continue;
            }
            if (!site_is_phased(matrix.sites[static_cast<size_t>(obs.site_index)])) continue;
            if (first_site < 0) first_site = obs.site_index;
            else union_parent(parent, first_site, obs.site_index);
        }
    }

    std::unordered_map<int, hts_pos_t> component_phase_set;
    for (size_t i = 0; i < matrix.sites.size(); ++i) {
        PhaseSite& site = matrix.sites[i];
        if (!site_is_phased(site)) continue;
        const int root = find_parent(parent, static_cast<int>(i));
        auto it = component_phase_set.find(root);
        if (it == component_phase_set.end() || site.order_pos < it->second) {
            component_phase_set[root] = site.order_pos;
        }
    }

    for (size_t i = 0; i < matrix.sites.size(); ++i) {
        PhaseSite& site = matrix.sites[i];
        if (!site_is_phased(site)) continue;
        site.phase_set = component_phase_set[find_parent(parent, static_cast<int>(i))];
    }
}

int update_site_consensus(PhaseMatrix& matrix) {
    std::vector<std::array<std::vector<int>, 3>> hap_counts(matrix.sites.size());
    for (size_t i = 0; i < matrix.sites.size(); ++i) {
        const int n = matrix.sites[i].n_alleles;
        hap_counts[i][1].assign(static_cast<size_t>(n), 0);
        hap_counts[i][2].assign(static_cast<size_t>(n), 0);
    }

    for (const ReadPhaseProfile& read : matrix.reads) {
        if (read.hap == 0) continue;
        for (const PhaseObservation& obs : read.observations) {
            if (obs.allele < 0 ||
                obs.site_index < 0 ||
                static_cast<size_t>(obs.site_index) >= matrix.sites.size()) {
                continue;
            }
            hap_counts[static_cast<size_t>(obs.site_index)]
                      [static_cast<size_t>(read.hap)]
                      [static_cast<size_t>(obs.allele)]++;
        }
    }

    int changed = 0;
    for (size_t site_i = 0; site_i < matrix.sites.size(); ++site_i) {
        PhaseSite& site = matrix.sites[site_i];
        if (!site_is_informative(site)) continue;
        for (int hap = 1; hap <= 2; ++hap) {
            int best = -1;
            int best_count = 0;
            for (int allele = 0; allele < site.n_alleles; ++allele) {
                const int count = hap_counts[site_i][hap][static_cast<size_t>(allele)];
                if (count > best_count) {
                    best_count = count;
                    best = allele;
                }
            }
            if (best != -1 && best != site.hap_to_cons_allele[static_cast<size_t>(hap)]) {
                site.hap_to_cons_allele[static_cast<size_t>(hap)] = best;
                changed = 1;
            }
        }
    }
    return changed;
}

int site_depth(const PhaseSite& site) {
    int depth = 0;
    for (int count : site.allele_counts) depth += count;
    return depth;
}

void write_counts(std::ostream& out, const std::vector<int>& counts) {
    for (size_t a = 0; a < counts.size(); ++a) {
        if (a > 0) out << ',';
        out << counts[a];
    }
}

void write_observations(std::ostream& out,
                        const PhaseMatrix& matrix,
                        const ReadPhaseProfile& read) {
    for (size_t i = 0; i < read.observations.size(); ++i) {
        if (i > 0) out << ',';
        const PhaseObservation& obs = read.observations[i];
        out << matrix.sites[static_cast<size_t>(obs.site_index)].site_id
            << ':' << obs.allele;
    }
}

void write_int_list_or_dot(std::ostream& out, const std::vector<int>& values) {
    if (values.empty()) {
        out << '.';
        return;
    }
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) out << ',';
        out << values[i];
    }
}

} // namespace

void phase_matrix(PhaseMatrix& matrix, int max_iters) {
    init_site_consensus(matrix);

    for (int iter = 0; iter < max_iters; ++iter) {
        int changed_reads = 0;
        for (ReadPhaseProfile& read : matrix.reads) {
            const int old_hap = read.hap;
            read.hap = assign_read_hap(matrix, read);
            if (read.hap != old_hap) ++changed_reads;
        }
        const int changed_sites = update_site_consensus(matrix);
        if (changed_reads == 0 && changed_sites == 0) break;
    }

    assign_phase_sets_from_read_components(matrix);

    for (ReadPhaseProfile& read : matrix.reads) {
        read.phase_set = (read.hap == 0) ? -1 : first_read_phase_set(matrix, read);
    }
}

void stitch_phase_chunks(std::vector<PhaseChunk>& chunks) {
    for (size_t chunk_i = 1; chunk_i < chunks.size(); ++chunk_i) {
        PhaseChunk& pre = chunks[chunk_i - 1];
        PhaseChunk& cur = chunks[chunk_i];
        if (pre.chrom != cur.chrom) continue;

        std::unordered_map<std::string, const ReadPhaseProfile*> pre_reads;
        pre_reads.reserve(pre.matrix.reads.size());
        for (const ReadPhaseProfile& read : pre.matrix.reads) {
            if (read.hap == 0 || read.phase_set <= 0) continue;
            pre_reads.emplace(read.read_name, &read);
        }

        struct StitchEvidence {
            int flip_hap_score = 0;
            hts_pos_t target_phase_set = -1;
            bool saw_overlap = false;
        };
        std::unordered_map<hts_pos_t, StitchEvidence> evidence_by_cur_ps;
        for (const ReadPhaseProfile& cur_read : cur.matrix.reads) {
            if (cur_read.hap == 0 || cur_read.phase_set <= 0) continue;
            auto it = pre_reads.find(cur_read.read_name);
            if (it == pre_reads.end()) continue;
            const ReadPhaseProfile& pre_read = *it->second;
            StitchEvidence& evidence = evidence_by_cur_ps[cur_read.phase_set];
            if (pre_read.hap == cur_read.hap) --evidence.flip_hap_score;
            else ++evidence.flip_hap_score;
            if (evidence.target_phase_set < pre_read.phase_set) {
                evidence.target_phase_set = pre_read.phase_set;
            }
            evidence.saw_overlap = true;
        }

        for (const auto& item : evidence_by_cur_ps) {
            const hts_pos_t cur_phase_set = item.first;
            const StitchEvidence& evidence = item.second;
            if (evidence.flip_hap_score == 0) continue;
            if (!evidence.saw_overlap || evidence.target_phase_set <= 0) continue;

            const bool do_flip = evidence.flip_hap_score > 0;
            for (PhaseSite& site : cur.matrix.sites) {
                if (site.phase_set != cur_phase_set) continue;
                if (do_flip) std::swap(site.hap_to_cons_allele[1], site.hap_to_cons_allele[2]);
                site.phase_set = evidence.target_phase_set;
            }
            for (ReadPhaseProfile& read : cur.matrix.reads) {
                if (read.phase_set != cur_phase_set) continue;
                if (do_flip && read.hap != 0) read.hap = 3 - read.hap;
                read.phase_set = evidence.target_phase_set;
            }
        }
    }
}

void write_phase_site_counts_tsv(std::ostream& out, const PhaseMatrix& matrix) {
    out << "SITE_INDEX\tSITE_ID\tCHROM\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\n";
    for (size_t i = 0; i < matrix.sites.size(); ++i) {
        const PhaseSite& site = matrix.sites[i];
        out << i << '\t'
            << site.site_id << '\t'
            << site.chrom << '\t'
            << site.order_pos << '\t'
            << site.n_alleles << '\t'
            << site_depth(site) << '\t';
        write_counts(out, site.allele_counts);
        out << '\n';
    }
}

void write_phase_read_profiles_tsv(std::ostream& out, const PhaseMatrix& matrix) {
    out << "READ\tN_OBS\tOBS\n";
    for (const ReadPhaseProfile& read : matrix.reads) {
        out << read.read_name << '\t' << read.observations.size() << '\t';
        write_observations(out, matrix, read);
        out << '\n';
    }
}

void write_phase_sites_tsv(std::ostream& out, const PhaseMatrix& matrix) {
    out << "SITE_INDEX\tSITE_ID\tCHROM\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\tPHASE_SET\tHAP1_ALLELE\tHAP2_ALLELE\tLEVEL\tPARENT\tROOT\tCOND_PARENT_ALLELES\n";
    for (size_t i = 0; i < matrix.sites.size(); ++i) {
        const PhaseSite& site = matrix.sites[i];
        out << i << '\t'
            << site.site_id << '\t'
            << site.chrom << '\t'
            << site.order_pos << '\t'
            << site.n_alleles << '\t'
            << site_depth(site) << '\t';
        write_counts(out, site.allele_counts);
        out << '\t'
            << site.phase_set << '\t'
            << site.hap_to_cons_allele[1] << '\t'
            << site.hap_to_cons_allele[2] << '\t'
            << site.level << '\t'
            << (site.parent_site_id.empty() ? "." : site.parent_site_id) << '\t'
            << (site.root_site_id.empty() ? "." : site.root_site_id) << '\t';
        write_int_list_or_dot(out, site.conditional_parent_alleles);
        out << '\n';
    }
}

void write_phase_reads_tsv(std::ostream& out, const PhaseMatrix& matrix) {
    out << "READ\tHAP\tPHASE_SET\tN_OBS\tOBS\n";
    for (const ReadPhaseProfile& read : matrix.reads) {
        out << read.read_name << '\t'
            << read.hap << '\t'
            << read.phase_set << '\t'
            << read.observations.size() << '\t';
        write_observations(out, matrix, read);
        out << '\n';
    }
}

} // namespace pgphase_collect
