#include "graph_bam_adapter.hpp"

#include "collect_output.hpp"
#include "collect_phase.hpp"

#include <algorithm>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

namespace {

struct GraphProfileObservation {
    int site_index = -1;
    int allele = -1;
};

std::string site_key(const GraphSite& site, size_t /*index*/) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + site.ref;
}

// Assigns a site to the chunk whose half-open interval [beg, end) contains the
// site's start position.  Using start-position (not overlap) ensures each site
// lands in exactly one chunk even when a snarl spans a chunk boundary.
bool site_starts_in_interval(const GraphSite& site,
                              const std::string& contig,
                              hts_pos_t beg,
                              hts_pos_t end) {
    const std::string site_contig = site.ref_contig.empty() ? site.chrom : site.ref_contig;
    if (!contig.empty() && site_contig != contig) return false;
    const hts_pos_t site_beg0 = (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
    return site_beg0 >= beg && site_beg0 < end;
}

bool graph_site_has_released_walk_storage(const GraphSite& site) {
    if (!site.eligible || site.allele_walks.size() < 2) return false;
    return std::all_of(site.allele_walks.begin(), site.allele_walks.end(),
                       [](const GraphWalk& walk) { return walk.empty(); });
}

void add_graph_candidate(GraphBamChunkBuildResult& out,
                         const GraphSite& site,
                         const std::string& key,
                         const std::vector<int>& allele_counts,
                         int tid) {
    CandidateVariant candidate;
    candidate.key.tid = tid;
    candidate.key.pos = site.order_pos() > 0 ? site.order_pos() : 1;
    candidate.key.type = VariantType::Snp;
    candidate.key.ref_len = 1;
    candidate.key.alt = key;
    candidate.counts.n_uniq_alles = static_cast<int>(allele_counts.size());
    candidate.counts.alle_covs = allele_counts;
    candidate.counts.ref_cov = allele_counts.empty() ? 0 : allele_counts[0];
    candidate.counts.alt_cov = 0;
    for (size_t allele = 1; allele < allele_counts.size(); ++allele) {
        candidate.counts.alt_cov += allele_counts[allele];
    }
    candidate.counts.total_cov = candidate.counts.ref_cov + candidate.counts.alt_cov;
    candidate.counts.forward_ref = candidate.counts.ref_cov;
    candidate.counts.forward_alt = candidate.counts.alt_cov;
    candidate.counts.allele_fraction =
        candidate.counts.total_cov > 0
            ? static_cast<double>(candidate.counts.alt_cov) / candidate.counts.total_cov
            : 0.0;
    candidate.counts.category = VariantCategory::CleanHetSnp;
    candidate.counts.candvarcate_initial = VariantCategory::CleanHetSnp;
    candidate.lcd_var_i_to_cate = kCandCleanHetSnp;
    candidate.hap_to_cons_alle[0] = -1;
    candidate.hap_to_cons_alle[1] = -1;
    candidate.hap_to_cons_alle[2] = -1;
    candidate.phase_set = -1;
    out.chunk.candidates.push_back(std::move(candidate));
    out.site_ids.push_back(key);
}

void add_read_profile(GraphBamChunkBuildResult& out,
                      const std::string& read_name,
                      const std::vector<GraphProfileObservation>& observations,
                      int mapq) {
    if (observations.empty()) return;

    ReadRecord read;
    read.tid = 0;
    read.input_index = 0;
    read.qname = read_name;
    read.mapq = mapq;
    read.is_skipped = false;
    read.beg = out.chunk.candidates[static_cast<size_t>(observations.front().site_index)].key.sort_pos();
    read.end = out.chunk.candidates[static_cast<size_t>(observations.back().site_index)].key.sort_pos();
    if (read.end < read.beg) read.end = read.beg;

    ReadVariantProfile profile;
    profile.read_id = static_cast<int>(out.chunk.reads.size());
    profile.start_var_idx = observations.front().site_index;
    profile.end_var_idx = observations.back().site_index;
    profile.alleles.assign(static_cast<size_t>(profile.end_var_idx - profile.start_var_idx + 1), -1);
    profile.alt_qi.assign(profile.alleles.size(), 0);
    for (const GraphProfileObservation& obs : observations) {
        const int offset = obs.site_index - profile.start_var_idx;
        if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
        profile.alleles[static_cast<size_t>(offset)] = obs.allele;
    }

    out.chunk.reads.push_back(std::move(read));
    out.chunk.read_var_profile.push_back(std::move(profile));
}

void rebuild_read_var_cr(BamChunk& chunk) {
    cgranges_t* cr = cr_init();
    if (cr == nullptr) throw std::runtime_error("failed to allocate graph read profile cgranges");
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.start_var_idx < 0 || profile.end_var_idx < profile.start_var_idx) continue;
        cr_add(cr, "cr", profile.start_var_idx, profile.end_var_idx + 1, profile.read_id);
    }
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

std::unordered_map<std::string, int> read_index_by_name(const BamChunk& chunk) {
    std::unordered_map<std::string, int> out;
    out.reserve(chunk.reads.size());
    for (size_t i = 0; i < chunk.reads.size(); ++i) {
        out.emplace(chunk.reads[i].qname, static_cast<int>(i));
    }
    return out;
}

void write_observations(std::ostream& out,
                        const GraphBamChunkBuildResult& graph_chunk,
                        const ReadVariantProfile& profile) {
    bool first = true;
    for (int site_i = profile.start_var_idx; site_i <= profile.end_var_idx; ++site_i) {
        const int offset = site_i - profile.start_var_idx;
        if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
        const int allele = profile.alleles[static_cast<size_t>(offset)];
        if (allele < 0) continue;
        if (!first) out << ',';
        first = false;
        const std::string& site_id = graph_chunk.site_ids[static_cast<size_t>(site_i)];
        out << site_id << ':' << allele;
    }
}

struct PhaseReadOutputRow {
    std::string read_name;
    int chunk_id = -1;
    int hap = 0;
    hts_pos_t phase_set = -1;
    bool has_phased_assignment = false;
    int copies = 0;
    int best_assignment_obs = -1;
    int assignment_chunk_id = std::numeric_limits<int>::max();
    std::unordered_map<std::string, int> allele_by_site;
};

void merge_phase_read_observation(PhaseReadOutputRow& row,
                                  const std::string& site_id,
                                  int allele) {
    auto [it, inserted] = row.allele_by_site.emplace(site_id, allele);
    if (!inserted && it->second != allele) it->second = -1;
}

void merge_phase_read_assignment(PhaseReadOutputRow& row,
                                 int chunk_id,
                                 int hap,
                                 hts_pos_t phase_set,
                                 int n_obs) {
    if (row.copies == 0) {
        row.chunk_id = chunk_id;
    } else if (row.chunk_id != chunk_id) {
        row.chunk_id = -1;
    }
    ++row.copies;

    if ((hap != 1 && hap != 2) || phase_set < 0) return;
    if (!row.has_phased_assignment ||
        n_obs > row.best_assignment_obs ||
        (n_obs == row.best_assignment_obs && chunk_id < row.assignment_chunk_id)) {
        row.hap = hap;
        row.phase_set = phase_set;
        row.has_phased_assignment = true;
        row.best_assignment_obs = n_obs;
        row.assignment_chunk_id = chunk_id;
    }
}

} // namespace

GraphBamChunkBuildResult build_graph_bam_chunk(const GraphSiteCatalog& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& contig,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts) {
    // Build a local contig→tid map.  The graph pipeline processes one contig at a
    // time, so this map always has a single entry and every site resolves to tid=0.
    // The map keeps key.tid and region.tid in sync and makes multi-contig extension
    // straightforward without any BAM header dependency.
    std::unordered_map<std::string, int> contig_tid_map;
    auto contig_tid = [&](const std::string& name) -> int {
        auto [it, inserted] = contig_tid_map.emplace(name, static_cast<int>(contig_tid_map.size()));
        return it->second;
    };
    const int chunk_tid = contig.empty() ? 0 : contig_tid(contig);

    GraphBamChunkBuildResult out;
    out.chunk.region.chunk_id = chunk_id;
    out.chunk.region.reg_chunk_i = chunk_id;
    out.chunk.region.tid = chunk_tid;
    out.chunk.region.beg = beg + 1;
    out.chunk.region.end = end;
    out.chunk.ref_beg = beg + 1;
    out.chunk.ref_end = end;
    out.chunk.chunk_min_qual = 60;
    out.chunk.chunk_first_quar_qual = 60;
    out.chunk.chunk_median_qual = 60;
    out.chunk.chunk_third_quar_qual = 60;
    out.chunk.chunk_max_qual = 60;
    out.chunk.up_ovlp_read_i.assign(1, {});
    out.chunk.down_ovlp_read_i.assign(1, {});
    out.chunk.n_up_ovlp_skip_reads.assign(1, 0);
    out.chunk.n_down_ovlp_skip_reads.assign(1, 0);

    std::unordered_map<std::string, int> site_to_candidate;
    std::vector<std::vector<int>> allele_counts;
    std::vector<int> parent_candidate;
    std::vector<std::vector<int>> conditional_parent_alleles;
    for (size_t site_i = 0; site_i < catalog.sites.size(); ++site_i) {
        const GraphSite& site = catalog.sites[site_i];
        if (!site_starts_in_interval(site, contig, beg, end)) continue;
        const bool released_walk_storage = graph_site_has_released_walk_storage(site);
        if (!released_walk_storage && !graph_site_is_queryable(site)) continue;
        const std::string key = site_key(site, site_i);
        const int n_alleles = static_cast<int>(site.allele_walks.size());
        if (n_alleles < 2) continue;
        site_to_candidate.emplace(key, static_cast<int>(out.chunk.candidates.size()));
        allele_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        parent_candidate.push_back(-1);
        conditional_parent_alleles.push_back(site.conditional_parent_alleles);
        const std::string& site_contig = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        const int site_tid = site_contig.empty() ? chunk_tid : contig_tid(site_contig);
        add_graph_candidate(out, site, key, allele_counts.back(), site_tid);
    }

    for (size_t site_i = 0; site_i < catalog.sites.size(); ++site_i) {
        const GraphSite& site = catalog.sites[site_i];
        if (site.parent.empty()) continue;
        const std::string key = site_key(site, site_i);
        auto child_it = site_to_candidate.find(key);
        auto parent_it = site_to_candidate.find(site.parent);
        if (child_it == site_to_candidate.end() || parent_it == site_to_candidate.end()) continue;
        parent_candidate[static_cast<size_t>(child_it->second)] = parent_it->second;
    }

    std::unordered_map<std::string, int> read_max_mapq;
    for (const GraphReadAllele& row : rows) {
        auto [it, inserted] = read_max_mapq.emplace(row.read_name, row.mapq);
        if (!inserted && row.mapq > it->second) it->second = row.mapq;
    }

    std::unordered_map<std::string, std::unordered_map<int, int>> allele_by_read_site;
    for (const GraphReadAllele& row : rows) {
        auto it = site_to_candidate.find(row.site_id);
        if (it == site_to_candidate.end()) continue;
        const int site_i = it->second;
        CandidateVariant& candidate = out.chunk.candidates[static_cast<size_t>(site_i)];
        if (row.allele < 0 || row.allele >= candidate.counts.n_uniq_alles) continue;
        std::unordered_map<int, int>& by_site = allele_by_read_site[row.read_name];
        auto obs_it = by_site.find(site_i);
        if (obs_it == by_site.end()) {
            by_site.emplace(site_i, row.allele);
        } else if (obs_it->second != row.allele) {
            obs_it->second = -1;
        }
    }

    std::unordered_map<std::string, std::vector<GraphProfileObservation>> read_obs;
    for (const auto& read_item : allele_by_read_site) {
        for (const auto& site_item : read_item.second) {
            const int site_i = site_item.first;
            const int allele = site_item.second;
            if (allele < 0) continue;
            const std::vector<int>& conditional = conditional_parent_alleles[static_cast<size_t>(site_i)];
            const int parent_i = parent_candidate[static_cast<size_t>(site_i)];
            if (!conditional.empty()) {
                if (parent_i < 0) continue;
                auto parent_obs = read_item.second.find(parent_i);
                if (parent_obs == read_item.second.end()) continue;
                if (std::find(conditional.begin(), conditional.end(), parent_obs->second) ==
                    conditional.end()) {
                    continue;
                }
            }
            ++allele_counts[static_cast<size_t>(site_i)][static_cast<size_t>(allele)];
            read_obs[read_item.first].push_back(GraphProfileObservation{site_i, allele});
        }
    }

    for (size_t site_i = 0; site_i < out.chunk.candidates.size(); ++site_i) {
        CandidateVariant& candidate = out.chunk.candidates[site_i];
        const std::vector<int>& counts = allele_counts[site_i];
        candidate.counts.alle_covs = counts;
        candidate.counts.ref_cov = counts.empty() ? 0 : counts[0];
        candidate.counts.alt_cov = 0;
        for (size_t allele = 1; allele < counts.size(); ++allele) candidate.counts.alt_cov += counts[allele];
        candidate.counts.total_cov = candidate.counts.ref_cov + candidate.counts.alt_cov;
        candidate.counts.forward_ref = candidate.counts.ref_cov;
        candidate.counts.forward_alt = candidate.counts.alt_cov;
        candidate.counts.allele_fraction =
            candidate.counts.total_cov > 0
                ? static_cast<double>(candidate.counts.alt_cov) / candidate.counts.total_cov
                : 0.0;
    }

    // Candidate filtering — three phases mirroring BAM §13 het classification.
    //
    // Phase 1: drop individual alt alleles below min_alt_depth (noise walks).
    // Build allele_remap[site_i][old_allele] = new_allele (or -1 if dropped).
    // Recompute per-site counts from surviving alleles only.
    const size_t n_cands = out.chunk.candidates.size();
    std::vector<std::vector<int>> allele_remap(n_cands);
    for (size_t i = 0; i < n_cands; ++i) {
        CandidateVariant& cand = out.chunk.candidates[i];
        std::vector<int>& ac = allele_counts[i];
        allele_remap[i].assign(ac.size(), -1);
        allele_remap[i][0] = 0;  // ref walk always kept at index 0
        std::vector<int> new_ac = {ac[0]};
        int next_allele = 1;
        for (size_t a = 1; a < ac.size(); ++a) {
            if (ac[a] >= opts.min_alt_depth) {
                allele_remap[i][a] = next_allele++;
                new_ac.push_back(ac[a]);
            }
        }
        allele_counts[i] = new_ac;
        cand.counts.alle_covs = new_ac;
        cand.counts.n_uniq_alles = static_cast<int>(new_ac.size());
        cand.counts.ref_cov = new_ac[0];
        cand.counts.alt_cov = 0;
        for (size_t a = 1; a < new_ac.size(); ++a) cand.counts.alt_cov += new_ac[a];
        cand.counts.total_cov = cand.counts.ref_cov + cand.counts.alt_cov;
        cand.counts.forward_ref = cand.counts.ref_cov;
        cand.counts.forward_alt = cand.counts.alt_cov;
        cand.counts.allele_fraction = cand.counts.total_cov > 0
            ? static_cast<double>(cand.counts.alt_cov) / cand.counts.total_cov
            : 0.0;
    }

    // Remap conditional_parent_alleles through parent's allele_remap so child
    // gates stay consistent after allele drops at the parent site.
    for (size_t i = 0; i < n_cands; ++i) {
        if (conditional_parent_alleles[i].empty()) continue;
        const int par = parent_candidate[i];
        if (par < 0 || static_cast<size_t>(par) >= n_cands) continue;
        const std::vector<int>& par_remap = allele_remap[static_cast<size_t>(par)];
        std::vector<int> new_cpa;
        for (int a : conditional_parent_alleles[i]) {
            if (a >= 0 && static_cast<size_t>(a) < par_remap.size() && par_remap[a] >= 0)
                new_cpa.push_back(par_remap[a]);
        }
        conditional_parent_alleles[i] = std::move(new_cpa);
    }

    // Build inverse allele mapping: for each surviving new allele index, record the original
    // allele-walk index so callers can map back to GraphSite.alts[orig-1].
    std::vector<std::vector<int>> allele_orig_idx(n_cands);
    for (size_t i = 0; i < n_cands; ++i) {
        allele_orig_idx[i].push_back(0);  // new index 0 = ref = orig index 0
        for (size_t old_a = 1; old_a < allele_remap[i].size(); ++old_a) {
            if (allele_remap[i][old_a] >= 0) {
                allele_orig_idx[i].push_back(static_cast<int>(old_a));
            }
        }
    }

    // Phase 2: decompose each surviving site into biallelic (ref vs alt_i) pairs,
    // mirroring the BAM path where every candidate is a single ref/alt pair.
    // AF and depth filters are applied per pair; filtered pairs go to out.filtered_sites.
    // A read observing ref contributes allele 0 to every pair from that site;
    // a read observing alt_i contributes allele 1 to only the pair for alt_i.
    struct NewPairEntry { int new_idx; int old_alt_phase1; };
    std::vector<std::vector<NewPairEntry>> old_site_to_pairs(n_cands);

    std::vector<CandidateVariant> new_cands;
    std::vector<std::string>      new_ids;
    std::vector<std::vector<int>> new_allele_counts;
    std::vector<int>              new_par;
    std::vector<std::vector<int>> new_cpa;
    std::vector<std::vector<int>> new_orig_idx;
    // Most sites are biallelic (1 surviving alt); reserve n_cands as a lower bound.
    new_cands.reserve(n_cands);
    new_ids.reserve(n_cands);
    new_allele_counts.reserve(n_cands);
    new_par.reserve(n_cands);
    new_cpa.reserve(n_cands);
    new_orig_idx.reserve(n_cands);

    for (size_t i = 0; i < n_cands; ++i) {
        const std::vector<int>& ac = allele_counts[i];
        if (ac.size() < 2) {
            out.filtered_sites.push_back({out.site_ids[i], ac[0], 0, ac[0], 0.0, "ref_only"});
            continue;
        }
        const int ref_c = ac[0];
        for (size_t a = 1; a < ac.size(); ++a) {
            const int alt_c   = ac[a];
            const int total_c = ref_c + alt_c;
            const double af   = total_c > 0 ? static_cast<double>(alt_c) / total_c : 0.0;
            const int orig_alt = allele_orig_idx[i][a];

            const std::string pair_id = (ac.size() > 2)
                ? out.site_ids[i] + ":" + std::to_string(orig_alt)
                : out.site_ids[i];

            std::string drop_reason;
            if (total_c < opts.min_depth) {
                drop_reason = "low_depth";
            } else if (af > opts.max_af) {
                drop_reason = "high_af";
            } else if (af < opts.min_af) {
                drop_reason = "low_af";
            }
            if (!drop_reason.empty()) {
                out.filtered_sites.push_back({pair_id, ref_c, alt_c, total_c, af,
                                              std::move(drop_reason)});
                continue;
            }

            const int new_idx = static_cast<int>(new_cands.size());
            old_site_to_pairs[i].push_back({new_idx, static_cast<int>(a)});

            CandidateVariant pair_cand = out.chunk.candidates[i];
            pair_cand.counts.alle_covs       = {ref_c, alt_c};
            pair_cand.counts.n_uniq_alles    = 2;
            pair_cand.counts.ref_cov         = ref_c;
            pair_cand.counts.alt_cov         = alt_c;
            pair_cand.counts.total_cov       = total_c;
            pair_cand.counts.forward_ref     = ref_c;
            pair_cand.counts.forward_alt     = alt_c;
            pair_cand.counts.allele_fraction = af;
            new_cands.push_back(pair_cand);
            new_ids.push_back(pair_id);
            new_allele_counts.push_back({ref_c, alt_c});
            new_par.push_back(-1);
            new_cpa.push_back({});
            new_orig_idx.push_back({0, orig_alt});
        }
    }

    // Remap parent/conditional-parent-allele indices into the new biallelic pair space.
    // After Phase 1, conditional_parent_alleles[i] holds phase-1 parent allele indices.
    for (size_t i = 0; i < n_cands; ++i) {
        const int old_par = parent_candidate[i];
        if (old_par < 0 || conditional_parent_alleles[i].empty()) continue;
        const auto& parent_pairs = old_site_to_pairs[static_cast<size_t>(old_par)];
        for (const NewPairEntry& child_entry : old_site_to_pairs[i]) {
            int assigned_par = -1;
            std::vector<int> cpa_entry;
            for (int p_phase1 : conditional_parent_alleles[i]) {
                for (const NewPairEntry& pe : parent_pairs) {
                    if (pe.old_alt_phase1 == p_phase1) {
                        assigned_par = pe.new_idx;
                        cpa_entry.push_back(1);
                        break;
                    }
                }
                if (assigned_par >= 0) break;
            }
            new_par[static_cast<size_t>(child_entry.new_idx)]  = assigned_par;
            new_cpa[static_cast<size_t>(child_entry.new_idx)]  = std::move(cpa_entry);
        }
    }

    out.chunk.candidates       = std::move(new_cands);
    out.site_ids               = std::move(new_ids);
    allele_counts              = std::move(new_allele_counts);
    parent_candidate           = std::move(new_par);
    conditional_parent_alleles = std::move(new_cpa);
    allele_orig_idx            = std::move(new_orig_idx);

    // Phase 3: remap read observations to the new biallelic pair space.
    // Allele 0 (ref) fans out to all new pairs from its original site (allele 0 in each).
    // Allele j (alt) maps to allele 1 in the one pair that holds alt j.
    for (auto& [read_name, obs] : read_obs) {
        std::vector<GraphProfileObservation> new_obs;
        new_obs.reserve(obs.size());  // biallelic sites: exact; multiallelic: slight under-reserve
        for (const auto& o : obs) {
            if (o.site_index < 0 || static_cast<size_t>(o.site_index) >= n_cands) continue;
            const size_t old_si = static_cast<size_t>(o.site_index);
            const auto& pairs = old_site_to_pairs[old_si];
            if (pairs.empty()) continue;
            // Fast path: biallelic site (one surviving pair) — no loop needed.
            if (pairs.size() == 1) {
                if (o.allele == 0) {
                    new_obs.push_back({pairs[0].new_idx, 0});
                } else if (o.allele > 0 &&
                           static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                    if (allele_remap[old_si][static_cast<size_t>(o.allele)] == pairs[0].old_alt_phase1)
                        new_obs.push_back({pairs[0].new_idx, 1});
                }
                continue;
            }
            // General path: multiallelic site with multiple surviving pairs.
            if (o.allele == 0) {
                for (const NewPairEntry& pe : pairs)
                    new_obs.push_back({pe.new_idx, 0});
            } else if (o.allele > 0 &&
                       static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                const int phase1_alt = allele_remap[old_si][static_cast<size_t>(o.allele)];
                if (phase1_alt >= 0) {
                    for (const NewPairEntry& pe : pairs) {
                        if (pe.old_alt_phase1 == phase1_alt) {
                            new_obs.push_back({pe.new_idx, 1});
                            break;
                        }
                    }
                }
            }
        }
        obs = std::move(new_obs);
    }

    std::vector<std::string> read_obs_order;
    read_obs_order.reserve(read_obs.size());
    for (const auto& item : read_obs) read_obs_order.push_back(item.first);
    std::sort(read_obs_order.begin(), read_obs_order.end());

    for (const std::string& read_name : read_obs_order) {
        std::vector<GraphProfileObservation>& observations = read_obs[read_name];
        std::sort(observations.begin(), observations.end(),
                  [](const GraphProfileObservation& lhs, const GraphProfileObservation& rhs) {
                      if (lhs.site_index != rhs.site_index) return lhs.site_index < rhs.site_index;
                      return lhs.allele < rhs.allele;
                  });
        std::vector<GraphProfileObservation> dedup;
        dedup.reserve(observations.size());
        for (const GraphProfileObservation& obs : observations) {
            if (!dedup.empty() && dedup.back().site_index == obs.site_index) {
                if (dedup.back().allele != obs.allele) dedup.back().allele = -1;
                continue;
            }
            dedup.push_back(obs);
        }
        dedup.erase(std::remove_if(dedup.begin(), dedup.end(),
                                   [](const GraphProfileObservation& obs) { return obs.allele < 0; }),
                    dedup.end());
        const auto mapq_it = read_max_mapq.find(read_name);
        const int mapq = mapq_it != read_max_mapq.end() ? mapq_it->second : 255;
        add_read_profile(out, read_name, dedup, mapq);
    }

    out.chunk.haps.assign(out.chunk.reads.size(), 0);
    out.chunk.phase_sets.assign(out.chunk.reads.size(), -1);
    rebuild_read_var_cr(out.chunk);
    out.site_allele_orig_idx = std::move(allele_orig_idx);
    return out;
}

void write_graph_bam_filtered_sites_tsv(std::ostream& out,
                                        const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    out << "site_id\tref_cov\talt_cov\ttotal_cov\tallele_fraction\tfilter_reason\n";
    for (const auto& gc : graph_chunks) {
        for (const auto& fs : gc.filtered_sites) {
            out << fs.site_id << '\t' << fs.ref_cov << '\t' << fs.alt_cov << '\t'
                << fs.total_cov << '\t' << fs.allele_fraction << '\t'
                << fs.filter_reason << '\n';
        }
    }
}

void populate_graph_chunk_overlaps(std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    for (size_t i = 1; i < graph_chunks.size(); ++i) {
        BamChunk& pre = graph_chunks[i - 1].chunk;
        BamChunk& cur = graph_chunks[i].chunk;
        pre.down_ovlp_read_i.assign(1, {});
        cur.up_ovlp_read_i.assign(1, {});
        pre.n_down_ovlp_skip_reads.assign(1, 0);
        cur.n_up_ovlp_skip_reads.assign(1, 0);
        const std::unordered_map<std::string, int> pre_reads = read_index_by_name(pre);
        const std::unordered_map<std::string, int> cur_reads = read_index_by_name(cur);
        for (const auto& item : cur_reads) {
            auto pre_it = pre_reads.find(item.first);
            if (pre_it == pre_reads.end()) continue;
            pre.down_ovlp_read_i[0].push_back(pre_it->second);
            cur.up_ovlp_read_i[0].push_back(item.second);
        }
    }
}

void phase_graph_bam_chunks(std::vector<GraphBamChunkBuildResult>& graph_chunks,
                            const Options& opts) {
    for (GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        assign_hap_based_on_germline_het_vars_kmeans(graph_chunk.chunk, opts, kCandHetVarCate);
    }
    populate_graph_chunk_overlaps(graph_chunks);
    std::vector<BamChunk> chunks;
    chunks.reserve(graph_chunks.size());
    for (GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        chunks.push_back(std::move(graph_chunk.chunk));
    }
    stitch_chunk_haps(chunks, &opts, nullptr);
    for (size_t i = 0; i < chunks.size(); ++i) {
        graph_chunks[i].chunk = std::move(chunks[i]);
    }
}

void write_graph_bam_site_counts_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\n";
    for (const GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        for (size_t i = 0; i < graph_chunk.chunk.candidates.size(); ++i) {
            const CandidateVariant& candidate = graph_chunk.chunk.candidates[i];
            out << graph_chunk.chunk.region.chunk_id << '\t'
                << i << '\t'
                << graph_chunk.site_ids[i] << '\t'
                << candidate.key.pos << '\t'
                << candidate.counts.n_uniq_alles << '\t'
                << candidate.counts.total_cov << '\t';
            for (size_t allele = 0; allele < candidate.counts.alle_covs.size(); ++allele) {
                if (allele > 0) out << ',';
                out << candidate.counts.alle_covs[allele];
            }
            out << '\n';
        }
    }
}

void write_graph_bam_read_profiles_tsv(std::ostream& out,
                                       const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    out << "CHUNK_ID\tREAD\tN_OBS\tOBS\n";
    for (const GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        for (const ReadVariantProfile& profile : graph_chunk.chunk.read_var_profile) {
            const ReadRecord& read = graph_chunk.chunk.reads[static_cast<size_t>(profile.read_id)];
            int n_obs = 0;
            for (int allele : profile.alleles) {
                if (allele >= 0) ++n_obs;
            }
            out << graph_chunk.chunk.region.chunk_id << '\t'
                << read.qname << '\t'
                << n_obs << '\t';
            write_observations(out, graph_chunk, profile);
            out << '\n';
        }
    }
}

void write_graph_bam_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\tPHASE_SET\tHAP1_ALLELE\tHAP2_ALLELE\n";
    for (const GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        for (size_t i = 0; i < graph_chunk.chunk.candidates.size(); ++i) {
            const CandidateVariant& candidate = graph_chunk.chunk.candidates[i];
            out << graph_chunk.chunk.region.chunk_id << '\t'
                << i << '\t'
                << graph_chunk.site_ids[i] << '\t'
                << candidate.key.pos << '\t'
                << candidate.counts.n_uniq_alles << '\t'
                << candidate.counts.total_cov << '\t';
            for (size_t allele = 0; allele < candidate.counts.alle_covs.size(); ++allele) {
                if (allele > 0) out << ',';
                out << candidate.counts.alle_covs[allele];
            }
            out << '\t'
                << candidate.phase_set << '\t'
                << candidate.hap_to_cons_alle[1] << '\t'
                << candidate.hap_to_cons_alle[2] << '\n';
        }
    }
}

void write_graph_bam_variants_tsv(std::ostream& out,
                                  const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                  const GraphSiteCatalog& catalog) {
    // Build site_id → GraphSite lookup. Biallelic pair IDs ("id:orig_alt") strip the suffix.
    std::unordered_map<std::string, const GraphSite*> site_map;
    for (const GraphSite& site : catalog.sites) {
        if (!site.id.empty() && site.id != ".")
            site_map.emplace(site.id, &site);
    }

    out << "CHROM\tPOS\tTYPE\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tLOW_QUAL_COUNT"
           "\tFORWARD_REF\tREVERSE_REF\tFORWARD_ALT\tREVERSE_ALT"
           "\tAF\tCATEGORY\tINIT_CAT\tPHASE_SET\tHAP_ALT\tHAP_REF\n";

    for (const GraphBamChunkBuildResult& gc : graph_chunks) {
        for (size_t ci = 0; ci < gc.chunk.candidates.size(); ++ci) {
            const CandidateVariant& cand = gc.chunk.candidates[ci];
            if (ci >= gc.site_ids.size()) continue;

            const std::string& raw_id = gc.site_ids[ci];
            auto it = site_map.find(raw_id);
            if (it == site_map.end()) {
                const auto colon_pos = raw_id.rfind(':');
                if (colon_pos != std::string::npos)
                    it = site_map.find(raw_id.substr(0, colon_pos));
                if (it == site_map.end()) continue;
            }
            const GraphSite& site = *it->second;

            const std::string& chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
            const std::string& ref_seq = site.ref;
            if (ref_seq.empty()) continue;

            // Map the surviving alt allele back to the original VCF ALT sequence.
            int orig_walk_idx = 1;
            if (ci < gc.site_allele_orig_idx.size() &&
                gc.site_allele_orig_idx[ci].size() > 1) {
                orig_walk_idx = gc.site_allele_orig_idx[ci][1];
            }
            const int alt_seq_idx = orig_walk_idx - 1;
            if (alt_seq_idx < 0 || alt_seq_idx >= static_cast<int>(site.alts.size())) continue;
            const std::string& alt_seq = site.alts[static_cast<size_t>(alt_seq_idx)];
            if (alt_seq.empty() || alt_seq == "*") continue;

            VariantType vtype;
            if (alt_seq.size() > ref_seq.size()) {
                vtype = VariantType::Insertion;
            } else if (ref_seq.size() > alt_seq.size()) {
                vtype = VariantType::Deletion;
            } else {
                vtype = VariantType::Snp;
            }

            // Resolve hap assignments (-1 = unphased → treat as ref for output).
            int h1 = cand.hap_to_cons_alle[1];
            int h2 = cand.hap_to_cons_alle[2];
            if (h1 == -1 && h2 == -1) h1 = h2 = cand.hap_to_cons_alle[0];
            if (h1 < 0) h1 = 0;
            if (h2 < 0) h2 = 0;

            out << chrom << '\t' << site.pos << '\t' << type_name(vtype) << '\t'
                << ref_seq << '\t' << alt_seq << '\t'
                << cand.counts.total_cov << '\t' << cand.counts.ref_cov << '\t'
                << cand.counts.alt_cov << '\t' << cand.counts.low_qual_cov << '\t'
                << cand.counts.forward_ref << '\t' << cand.counts.reverse_ref << '\t'
                << cand.counts.forward_alt << '\t' << cand.counts.reverse_alt << '\t'
                << cand.counts.allele_fraction << '\t'
                << category_name(cand.counts.category) << '\t'
                << category_name(cand.counts.candvarcate_initial) << '\t'
                << cand.phase_set << '\t'
                << (h1 > 0 ? 1 : 0) << '\t' << (h2 > 0 ? 1 : 0) << '\n';
        }
    }
}

void write_graph_bam_phase_reads_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                     const std::vector<std::string>& all_read_names) {
    out << "CHUNK_ID\tREAD\tHAP\tPHASE_SET\tN_OBS\tOBS\n";

    std::unordered_map<std::string, PhaseReadOutputRow> rows_by_read;
    for (const GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        const BamChunk& chunk = graph_chunk.chunk;
        for (const ReadVariantProfile& profile : chunk.read_var_profile) {
            const size_t read_i = static_cast<size_t>(profile.read_id);
            const ReadRecord& read = chunk.reads[read_i];
            PhaseReadOutputRow& row = rows_by_read[read.qname];
            if (row.read_name.empty()) row.read_name = read.qname;

            const int hap = read_i < chunk.haps.size() ? chunk.haps[read_i] : 0;
            const hts_pos_t phase_set =
                read_i < chunk.phase_sets.size()
                    ? chunk.phase_sets[read_i]
                    : static_cast<hts_pos_t>(-1);

            int n_obs = 0;
            for (int site_i = profile.start_var_idx; site_i <= profile.end_var_idx; ++site_i) {
                const int offset = site_i - profile.start_var_idx;
                if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
                const int allele = profile.alleles[static_cast<size_t>(offset)];
                if (allele < 0) continue;
                ++n_obs;
                if (site_i < 0 || static_cast<size_t>(site_i) >= graph_chunk.site_ids.size()) continue;
                merge_phase_read_observation(row,
                                             graph_chunk.site_ids[static_cast<size_t>(site_i)],
                                             allele);
            }
            merge_phase_read_assignment(row, chunk.region.chunk_id, hap, phase_set, n_obs);
        }
    }

    for (const std::string& name : all_read_names) {
        auto [it, inserted] = rows_by_read.emplace(name, PhaseReadOutputRow{});
        if (inserted) it->second.read_name = name;
    }

    std::vector<std::string> read_names;
    read_names.reserve(rows_by_read.size());
    for (const auto& item : rows_by_read) read_names.push_back(item.first);
    std::sort(read_names.begin(), read_names.end());

    for (const std::string& read_name : read_names) {
        PhaseReadOutputRow& row = rows_by_read[read_name];
        if (!row.has_phased_assignment) {
            row.hap = 0;
            row.phase_set = -1;
        }

        std::vector<std::string> obs_keys;
        obs_keys.reserve(row.allele_by_site.size());
        for (const auto& obs : row.allele_by_site) {
            if (obs.second >= 0) obs_keys.push_back(obs.first);
        }
        std::sort(obs_keys.begin(), obs_keys.end());

        out << row.chunk_id << '\t'
            << read_name << '\t'
            << row.hap << '\t'
            << row.phase_set << '\t'
            << obs_keys.size() << '\t';
        for (size_t i = 0; i < obs_keys.size(); ++i) {
            if (i > 0) out << ',';
            const std::string& site_id = obs_keys[i];
            out << site_id << ':' << row.allele_by_site[site_id];
        }
        out << '\n';
    }
}

} // namespace pgphase_collect
