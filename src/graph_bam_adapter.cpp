#include "graph_bam_adapter.hpp"

#include "collect_output.hpp"
#include "collect_phase.hpp"
#include "fisher_exact.hpp"

#include <algorithm>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern "C" {
#include "cgranges.h"
#include "htslib/sam.h"
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

void merge_phase_read_observation(PhaseReadOutputRow& row,
                                  const std::string& site_id,
                                  int allele) {
    auto [it, inserted] = row.allele_by_site.emplace(site_id, allele);
    if (!inserted && it->second != allele) it->second = -1;
}

void merge_phase_read_assignment(PhaseReadOutputRow& row,
                                 int chunk_id,
                                 int hap,
                                 hts_pos_t phase_set) {
    if (row.copies == 0) {
        row.chunk_id = chunk_id;
    } else if (row.chunk_id != chunk_id) {
        row.chunk_id = -1;
    }
    ++row.copies;

    // Match collect_bam_output PhasedAlignmentWriter::write_chunks: HP/PS come from the one chunk
    // that owns the read at emission time. Overlap reads skip upstream chunks there; here we merge
    // finalized chunks in order, so the latest visit matches downstream ownership / stitched state.
    row.hap                   = hap;
    row.phase_set             = phase_set;
    row.has_phased_assignment = ((hap == 1 || hap == 2) && phase_set >= 0);
}

// Classify graph candidates using count-based rules matching the BAM pipeline's
// classify_variant_initial (collect_var.cpp), minus reference-sequence-dependent
// checks (homopolymer/repeat detection).  Runs in a single pass over candidates.
void classify_graph_candidates(BamChunk& chunk, const Options& opts) {
    for (auto& cand : chunk.candidates) {
        VariantCounts& c = cand.counts;

        // Low coverage / low alt depth.
        if (c.total_cov < opts.min_depth || c.alt_cov < opts.min_alt_depth) {
            c.category = VariantCategory::LowCoverage;
            c.candvarcate_initial = VariantCategory::LowCoverage;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::LowCoverage);
            continue;
        }

        // ONT strand-bias filter (Fisher exact test on alt forward vs reverse).
        if (opts.is_ont()) {
            const int fa = c.forward_alt;
            const int ra = c.reverse_alt;
            const int expected = (fa + ra) / 2;
            if (expected > 0) {
                const double p = fisher_exact_two_tail(fa, ra, expected, expected);
                if (p < opts.strand_bias_pval) {
                    c.category = VariantCategory::StrandBias;
                    c.candvarcate_initial = VariantCategory::StrandBias;
                    cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::StrandBias);
                    continue;
                }
            }
        }

        // Low allele fraction → folded to LowCoverage (matches BAM pipeline pass 2).
        if (c.allele_fraction < opts.min_af) {
            c.category = VariantCategory::LowCoverage;
            c.candvarcate_initial = VariantCategory::LowAlleleFraction;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::LowCoverage);
            continue;
        }

        // Homozygous (high AF).
        if (c.allele_fraction > opts.max_af) {
            c.category = VariantCategory::CleanHom;
            c.candvarcate_initial = VariantCategory::CleanHom;
            cand.lcd_var_i_to_cate = category_to_flag(VariantCategory::CleanHom);
            continue;
        }

        // Surviving het — type-aware classification.
        if (cand.key.type == VariantType::Snp) {
            c.category = VariantCategory::CleanHetSnp;
            c.candvarcate_initial = VariantCategory::CleanHetSnp;
            cand.lcd_var_i_to_cate = kCandCleanHetSnp;
        } else {
            c.category = VariantCategory::CleanHetIndel;
            c.candvarcate_initial = VariantCategory::CleanHetIndel;
            cand.lcd_var_i_to_cate = kCandCleanHetIndel;
        }
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
    out.graph_phase_contig = contig;
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
    // Per-site, per-allele strand counts — accumulated from raw rows before dedup.
    std::vector<std::vector<int>> fwd_strand_counts;
    std::vector<std::vector<int>> rev_strand_counts;
    std::vector<int> parent_candidate;
    std::vector<std::vector<int>> conditional_parent_alleles;
    std::vector<size_t> catalog_site_idx_phase1;
    for (size_t site_i = 0; site_i < catalog.sites.size(); ++site_i) {
        const GraphSite& site = catalog.sites[site_i];
        if (!site_starts_in_interval(site, contig, beg, end)) continue;
        const std::string sid = site_key(site, site_i);
        if (!site.eligible) {
            out.filtered_sites.push_back({sid, 0, 0, 0, 0.0, "precandidate_ineligible"});
            continue;
        }
        const bool released_walk_storage = graph_site_has_released_walk_storage(site);
        if (!released_walk_storage && !graph_site_is_queryable(site)) {
            out.filtered_sites.push_back({sid, 0, 0, 0, 0.0, "precandidate_not_queryable"});
            continue;
        }
        const int n_alleles = static_cast<int>(site.allele_walks.size());
        if (n_alleles < 2) {
            out.filtered_sites.push_back({sid, 0, 0, 0, 0.0, "precandidate_monoallelic"});
            continue;
        }
        site_to_candidate.emplace(sid, static_cast<int>(out.chunk.candidates.size()));
        allele_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        fwd_strand_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        rev_strand_counts.emplace_back(static_cast<size_t>(n_alleles), 0);
        parent_candidate.push_back(-1);
        conditional_parent_alleles.push_back(site.conditional_parent_alleles);
        const std::string& site_contig = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        const int site_tid = site_contig.empty() ? chunk_tid : contig_tid(site_contig);
        catalog_site_idx_phase1.push_back(site_i);
        add_graph_candidate(out, site, sid, allele_counts.back(), site_tid);
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

    // Pack allele + reverse into a single int: low 16 bits = allele, bit 16 = reverse.
    // Conflict sentinel is -1 (allele < 0).
    constexpr int kRevBit = 0x10000;
    auto pack_allele_rev = [](int allele, bool rev) -> int {
        return allele | (rev ? kRevBit : 0);
    };
    auto unpack_allele = [](int packed) -> int { return packed & 0xFFFF; };
    auto unpack_reverse = [&](int packed) -> bool { return (packed & kRevBit) != 0; };

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
            by_site.emplace(site_i, pack_allele_rev(row.allele, row.reverse));
        } else if (obs_it->second >= 0 && unpack_allele(obs_it->second) != row.allele) {
            obs_it->second = -1;  // conflict
        }
    }

    std::unordered_map<std::string, std::vector<GraphProfileObservation>> read_obs;
    for (const auto& read_item : allele_by_read_site) {
        for (const auto& site_item : read_item.second) {
            const int site_i = site_item.first;
            const int packed = site_item.second;
            if (packed < 0) continue;
            const int allele = unpack_allele(packed);
            const bool rev = unpack_reverse(packed);
            const std::vector<int>& conditional = conditional_parent_alleles[static_cast<size_t>(site_i)];
            const int parent_i = parent_candidate[static_cast<size_t>(site_i)];
            if (!conditional.empty()) {
                if (parent_i < 0) continue;
                auto parent_obs = read_item.second.find(parent_i);
                if (parent_obs == read_item.second.end() || parent_obs->second < 0) continue;
                const int parent_allele = unpack_allele(parent_obs->second);
                if (std::find(conditional.begin(), conditional.end(), parent_allele) ==
                    conditional.end()) {
                    continue;
                }
            }
            ++allele_counts[static_cast<size_t>(site_i)][static_cast<size_t>(allele)];
            auto& sc = rev ? rev_strand_counts[static_cast<size_t>(site_i)]
                           : fwd_strand_counts[static_cast<size_t>(site_i)];
            ++sc[static_cast<size_t>(allele)];
            read_obs[read_item.first].push_back(GraphProfileObservation{site_i, allele});
        }
    }

    for (size_t site_i = 0; site_i < out.chunk.candidates.size(); ++site_i) {
        CandidateVariant& candidate = out.chunk.candidates[site_i];
        const std::vector<int>& counts = allele_counts[site_i];
        const std::vector<int>& fwd = fwd_strand_counts[site_i];
        const std::vector<int>& rev = rev_strand_counts[site_i];
        candidate.counts.alle_covs = counts;
        candidate.counts.ref_cov = counts.empty() ? 0 : counts[0];
        candidate.counts.alt_cov = 0;
        for (size_t allele = 1; allele < counts.size(); ++allele) candidate.counts.alt_cov += counts[allele];
        candidate.counts.total_cov = candidate.counts.ref_cov + candidate.counts.alt_cov;
        candidate.counts.forward_ref = fwd.empty() ? 0 : fwd[0];
        candidate.counts.reverse_ref = rev.empty() ? 0 : rev[0];
        candidate.counts.forward_alt = 0;
        candidate.counts.reverse_alt = 0;
        for (size_t allele = 1; allele < fwd.size(); ++allele) candidate.counts.forward_alt += fwd[allele];
        for (size_t allele = 1; allele < rev.size(); ++allele) candidate.counts.reverse_alt += rev[allele];
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
        std::vector<int>& fc = fwd_strand_counts[i];
        std::vector<int>& rc = rev_strand_counts[i];
        allele_remap[i].assign(ac.size(), -1);
        allele_remap[i][0] = 0;  // ref walk always kept at index 0
        std::vector<int> new_ac = {ac[0]};
        std::vector<int> new_fc = {fc[0]};
        std::vector<int> new_rc = {rc[0]};
        int next_allele = 1;
        for (size_t a = 1; a < ac.size(); ++a) {
            if (ac[a] >= opts.min_alt_depth) {
                allele_remap[i][a] = next_allele++;
                new_ac.push_back(ac[a]);
                new_fc.push_back(a < fc.size() ? fc[a] : 0);
                new_rc.push_back(a < rc.size() ? rc[a] : 0);
            }
        }
        allele_counts[i] = new_ac;
        fwd_strand_counts[i] = new_fc;
        rev_strand_counts[i] = new_rc;
        cand.counts.alle_covs = new_ac;
        cand.counts.n_uniq_alles = static_cast<int>(new_ac.size());
        cand.counts.ref_cov = new_ac[0];
        cand.counts.alt_cov = 0;
        for (size_t a = 1; a < new_ac.size(); ++a) cand.counts.alt_cov += new_ac[a];
        cand.counts.total_cov = cand.counts.ref_cov + cand.counts.alt_cov;
        cand.counts.forward_ref = new_fc[0];
        cand.counts.reverse_ref = new_rc[0];
        cand.counts.forward_alt = 0;
        cand.counts.reverse_alt = 0;
        for (size_t a = 1; a < new_fc.size(); ++a) cand.counts.forward_alt += new_fc[a];
        for (size_t a = 1; a < new_rc.size(); ++a) cand.counts.reverse_alt += new_rc[a];
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
    std::vector<std::vector<int>> new_fwd_strand;
    std::vector<std::vector<int>> new_rev_strand;
    std::vector<int>              new_par;
    std::vector<std::vector<int>> new_cpa;
    std::vector<std::vector<int>> new_orig_idx;
    // Most sites are biallelic (1 surviving alt); reserve n_cands as a lower bound.
    new_cands.reserve(n_cands);
    new_ids.reserve(n_cands);
    new_allele_counts.reserve(n_cands);
    new_fwd_strand.reserve(n_cands);
    new_rev_strand.reserve(n_cands);
    new_par.reserve(n_cands);
    new_cpa.reserve(n_cands);
    new_orig_idx.reserve(n_cands);

    for (size_t i = 0; i < n_cands; ++i) {
        const std::vector<int>& ac = allele_counts[i];
        if (ac.size() < 2) {
            out.filtered_sites.push_back({out.site_ids[i], ac[0], 0, ac[0], 0.0, "ref_only"});
            continue;
        }
        const std::vector<int>& fc = fwd_strand_counts[i];
        const std::vector<int>& rc_s = rev_strand_counts[i];
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

            const int fwd_ref = fc.empty() ? 0 : fc[0];
            const int rev_ref = rc_s.empty() ? 0 : rc_s[0];
            const int fwd_alt = a < fc.size() ? fc[a] : 0;
            const int rev_alt = a < rc_s.size() ? rc_s[a] : 0;

            const int new_idx = static_cast<int>(new_cands.size());
            old_site_to_pairs[i].push_back({new_idx, static_cast<int>(a)});

            CandidateVariant pair_cand = out.chunk.candidates[i];

            // Derive variant type from VCF REF/ALT sequence lengths.
            {
                const GraphSite& site = catalog.sites[catalog_site_idx_phase1[i]];
                const size_t ref_len = site.ref.size();
                const size_t alt_idx = static_cast<size_t>(orig_alt - 1);
                const size_t alt_len = alt_idx < site.alts.size() ? site.alts[alt_idx].size() : ref_len;
                if (ref_len < alt_len) {
                    pair_cand.key.type = VariantType::Insertion;
                    pair_cand.key.ref_len = 1;
                } else if (ref_len > alt_len) {
                    pair_cand.key.type = VariantType::Deletion;
                    pair_cand.key.ref_len = static_cast<int>(ref_len);
                } else {
                    pair_cand.key.type = VariantType::Snp;
                    pair_cand.key.ref_len = static_cast<int>(ref_len);
                }
            }

            pair_cand.counts.alle_covs       = {ref_c, alt_c};
            pair_cand.counts.n_uniq_alles    = 2;
            pair_cand.counts.ref_cov         = ref_c;
            pair_cand.counts.alt_cov         = alt_c;
            pair_cand.counts.total_cov       = total_c;
            pair_cand.counts.forward_ref     = fwd_ref;
            pair_cand.counts.reverse_ref     = rev_ref;
            pair_cand.counts.forward_alt     = fwd_alt;
            pair_cand.counts.reverse_alt     = rev_alt;
            pair_cand.counts.allele_fraction = af;
            new_cands.push_back(pair_cand);
            new_ids.push_back(pair_id);
            new_allele_counts.push_back({ref_c, alt_c});
            new_fwd_strand.push_back({fwd_ref, fwd_alt});
            new_rev_strand.push_back({rev_ref, rev_alt});
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
    fwd_strand_counts          = std::move(new_fwd_strand);
    rev_strand_counts          = std::move(new_rev_strand);
    parent_candidate           = std::move(new_par);
    conditional_parent_alleles = std::move(new_cpa);
    allele_orig_idx            = std::move(new_orig_idx);

    // Classify candidates before building read profiles so that pruned sites
    // (LowCoverage, StrandBias) never enter the profile / cr_overlap index.
    classify_graph_candidates(out.chunk, opts);

    // Build a fast lookup for pruned candidates.
    const size_t n_final_cands = out.chunk.candidates.size();
    std::vector<bool> cand_pruned(n_final_cands, false);
    for (size_t ci = 0; ci < n_final_cands; ++ci) {
        const VariantCategory cat = out.chunk.candidates[ci].counts.category;
        if (cat == VariantCategory::LowCoverage || cat == VariantCategory::StrandBias ||
            cat == VariantCategory::NonVariant) {
            cand_pruned[ci] = true;
        }
    }

    // Phase 3: remap read observations to the new biallelic pair space.
    // Allele 0 (ref) fans out to all new pairs from its original site (allele 0 in each).
    // Allele j (alt) maps to allele 1 in the one pair that holds alt j.
    // Observations at pruned candidates are dropped.
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
                if (cand_pruned[static_cast<size_t>(pairs[0].new_idx)]) continue;
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
                    if (!cand_pruned[static_cast<size_t>(pe.new_idx)])
                        new_obs.push_back({pe.new_idx, 0});
            } else if (o.allele > 0 &&
                       static_cast<size_t>(o.allele) < allele_remap[old_si].size()) {
                const int phase1_alt = allele_remap[old_si][static_cast<size_t>(o.allele)];
                if (phase1_alt >= 0) {
                    for (const NewPairEntry& pe : pairs) {
                        if (pe.old_alt_phase1 == phase1_alt) {
                            if (!cand_pruned[static_cast<size_t>(pe.new_idx)])
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

    constexpr size_t kNoPhase1 = std::numeric_limits<size_t>::max();
    std::vector<size_t> final_to_phase1(out.chunk.candidates.size(), kNoPhase1);
    for (size_t old_i = 0; old_i < old_site_to_pairs.size(); ++old_i) {
        for (const auto& pe : old_site_to_pairs[old_i]) {
            if (pe.new_idx >= 0 &&
                static_cast<size_t>(pe.new_idx) < final_to_phase1.size())
                final_to_phase1[static_cast<size_t>(pe.new_idx)] = old_i;
        }
    }

    out.variant_emit_rows.assign(out.chunk.candidates.size(), {});
    for (size_t ci = 0; ci < out.chunk.candidates.size(); ++ci) {
        const size_t pi = final_to_phase1[ci];
        if (pi == kNoPhase1 || pi >= catalog_site_idx_phase1.size()) continue;
        const GraphSite& site = catalog.sites[catalog_site_idx_phase1[pi]];
        GraphVariantEmitRow row;
        row.chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        row.pos = site.pos;
        row.snarl_id = ci < out.site_ids.size() ? out.site_ids[ci] : ".";
        row.ref = site.ref;

        int orig_walk_idx = 1;
        if (ci < out.site_allele_orig_idx.size() &&
            out.site_allele_orig_idx[ci].size() > 1) {
            orig_walk_idx = out.site_allele_orig_idx[ci][1];
        }
        const int alt_seq_idx = orig_walk_idx - 1;
        if (alt_seq_idx >= 0 && static_cast<size_t>(alt_seq_idx) < site.alts.size())
            row.alt = site.alts[static_cast<size_t>(alt_seq_idx)];
        out.variant_emit_rows[ci] = std::move(row);
    }

    return out;
}

void write_graph_bam_filtered_sites_tsv_header(std::ostream& out) {
    out << "site_id\tref_cov\talt_cov\ttotal_cov\tallele_fraction\tfilter_reason\n";
}

void write_graph_bam_filtered_sites_tsv_rows(std::ostream& out,
                                             const GraphBamChunkBuildResult& gc) {
    for (const auto& fs : gc.filtered_sites) {
        out << fs.site_id << '\t' << fs.ref_cov << '\t' << fs.alt_cov << '\t'
            << fs.total_cov << '\t' << fs.allele_fraction << '\t'
            << fs.filter_reason << '\n';
    }
}

void write_graph_bam_filtered_sites_tsv(std::ostream& out,
                                        const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    write_graph_bam_filtered_sites_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_bam_filtered_sites_tsv_rows(out, gc);
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

void populate_graph_chunk_pair_overlaps(GraphBamChunkBuildResult& pre,
                                        GraphBamChunkBuildResult& cur) {
    BamChunk& pre_c = pre.chunk;
    BamChunk& cur_c = cur.chunk;
    pre_c.down_ovlp_read_i.assign(1, {});
    cur_c.up_ovlp_read_i.assign(1, {});
    pre_c.n_down_ovlp_skip_reads.assign(1, 0);
    cur_c.n_up_ovlp_skip_reads.assign(1, 0);
    const std::unordered_map<std::string, int> pre_reads = read_index_by_name(pre_c);
    const std::unordered_map<std::string, int> cur_reads = read_index_by_name(cur_c);
    for (const auto& item : cur_reads) {
        auto pre_it = pre_reads.find(item.first);
        if (pre_it == pre_reads.end()) continue;
        pre_c.down_ovlp_read_i[0].push_back(pre_it->second);
        cur_c.up_ovlp_read_i[0].push_back(item.second);
    }
}

void assign_graph_chunk_hap(GraphBamChunkBuildResult& gc, const Options& opts) {
    assign_hap_based_on_germline_het_vars_kmeans(gc.chunk, opts, kCandGermlineClean);
}

void stitch_graph_chunk_pair(GraphBamChunkBuildResult& pre,
                             GraphBamChunkBuildResult& cur,
                             const Options& opts) {
    populate_graph_chunk_pair_overlaps(pre, cur);
    std::vector<BamChunk> pair;
    pair.reserve(2);
    pair.push_back(std::move(pre.chunk));
    pair.push_back(std::move(cur.chunk));
    stitch_chunk_haps(pair, &opts, nullptr);
    pre.chunk = std::move(pair[0]);
    cur.chunk = std::move(pair[1]);
}

void phase_graph_bam_chunks(std::vector<GraphBamChunkBuildResult>& graph_chunks,
                            const Options& opts) {
    for (GraphBamChunkBuildResult& graph_chunk : graph_chunks) {
        assign_hap_based_on_germline_het_vars_kmeans(graph_chunk.chunk, opts, kCandGermlineClean);
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

void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphBamChunkBuildResult& gc) {
    const BamChunk& chunk = gc.chunk;
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

        for (int site_i = profile.start_var_idx; site_i <= profile.end_var_idx; ++site_i) {
            const int offset = site_i - profile.start_var_idx;
            if (offset < 0 || static_cast<size_t>(offset) >= profile.alleles.size()) continue;
            const int allele = profile.alleles[static_cast<size_t>(offset)];
            if (allele < 0) continue;
            if (site_i < 0 || static_cast<size_t>(site_i) >= gc.site_ids.size()) continue;
            merge_phase_read_observation(row,
                                         gc.site_ids[static_cast<size_t>(site_i)],
                                         allele);
        }
        merge_phase_read_assignment(row, chunk.region.chunk_id, hap, phase_set);
    }
}

void write_graph_bam_site_counts_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\n";
}

void write_graph_bam_site_counts_tsv_rows(std::ostream& out,
                                          const GraphBamChunkBuildResult& gc) {
    for (size_t i = 0; i < gc.chunk.candidates.size(); ++i) {
        const CandidateVariant& candidate = gc.chunk.candidates[i];
        out << gc.chunk.region.chunk_id << '\t'
            << i << '\t'
            << gc.site_ids[i] << '\t'
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

void write_graph_bam_site_counts_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    write_graph_bam_site_counts_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_bam_site_counts_tsv_rows(out, gc);
}

void write_graph_bam_read_profiles_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tREAD\tN_OBS\tOBS\n";
}

void write_graph_bam_read_profiles_tsv_rows(std::ostream& out,
                                            const GraphBamChunkBuildResult& gc) {
    for (const ReadVariantProfile& profile : gc.chunk.read_var_profile) {
        const ReadRecord& read = gc.chunk.reads[static_cast<size_t>(profile.read_id)];
        int n_obs = 0;
        for (int allele : profile.alleles) {
            if (allele >= 0) ++n_obs;
        }
        out << gc.chunk.region.chunk_id << '\t'
            << read.qname << '\t'
            << n_obs << '\t';
        write_observations(out, gc, profile);
        out << '\n';
    }
}

void write_graph_bam_read_profiles_tsv(std::ostream& out,
                                       const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    write_graph_bam_read_profiles_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_bam_read_profiles_tsv_rows(out, gc);
}

void write_graph_bam_phase_sites_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\tPHASE_SET\tHAP1_ALLELE\tHAP2_ALLELE\n";
}

void write_graph_bam_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphBamChunkBuildResult& gc) {
    for (size_t i = 0; i < gc.chunk.candidates.size(); ++i) {
        const CandidateVariant& candidate = gc.chunk.candidates[i];
        out << gc.chunk.region.chunk_id << '\t'
            << i << '\t'
            << gc.site_ids[i] << '\t'
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

void write_graph_bam_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphBamChunkBuildResult>& graph_chunks) {
    write_graph_bam_phase_sites_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_bam_phase_sites_tsv_rows(out, gc);
}

void write_graph_bam_variants_tsv_header(std::ostream& out) {
    out << "CHROM\tPOS\tSNARL_ID\tTYPE\tREF\tALT\tDP\tREF_COUNT\tALT_COUNT\tLOW_QUAL_COUNT"
           "\tFORWARD_REF\tREVERSE_REF\tFORWARD_ALT\tREVERSE_ALT"
           "\tAF\tCATEGORY\tINIT_CAT\tPHASE_SET\tHAP_ALT\tHAP_REF\n";
}

static void write_graph_bam_variant_tsv_line(std::ostream& out,
                                             const CandidateVariant& cand,
                                             const std::string& chrom,
                                             hts_pos_t pos,
                                             const std::string& snarl_id,
                                             const std::string& ref_seq,
                                             const std::string& alt_seq) {
    if (ref_seq.empty()) return;
    if (alt_seq.empty() || alt_seq == "*") return;

    VariantType vtype;
    if (alt_seq.size() > ref_seq.size()) {
        vtype = VariantType::Insertion;
    } else if (ref_seq.size() > alt_seq.size()) {
        vtype = VariantType::Deletion;
    } else {
        vtype = VariantType::Snp;
    }

    int h1 = cand.hap_to_cons_alle[1];
    int h2 = cand.hap_to_cons_alle[2];
    if (h1 == -1 && h2 == -1) h1 = h2 = cand.hap_to_cons_alle[0];
    if (h1 < 0) h1 = 0;
    if (h2 < 0) h2 = 0;

    out << chrom << '\t' << pos << '\t' << snarl_id << '\t' << type_name(vtype) << '\t'
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

static void write_graph_bam_variant_one_row(std::ostream& out,
                                            const CandidateVariant& cand,
                                            const GraphSite& site,
                                            const GraphBamChunkBuildResult& gc,
                                            size_t ci) {
    const std::string& chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
    const std::string& ref_seq = site.ref;
    if (ref_seq.empty()) return;

    int orig_walk_idx = 1;
    if (ci < gc.site_allele_orig_idx.size() &&
        gc.site_allele_orig_idx[ci].size() > 1) {
        orig_walk_idx = gc.site_allele_orig_idx[ci][1];
    }
    const int alt_seq_idx = orig_walk_idx - 1;
    if (alt_seq_idx < 0 || alt_seq_idx >= static_cast<int>(site.alts.size())) return;
    const std::string& alt_seq = site.alts[static_cast<size_t>(alt_seq_idx)];
    const std::string& snarl_id = ci < gc.site_ids.size() ? gc.site_ids[ci] : std::string(".");
    write_graph_bam_variant_tsv_line(out, cand, chrom, site.pos, snarl_id, ref_seq, alt_seq);
}

void write_graph_bam_variants_tsv_rows(
    std::ostream& out,
    const GraphBamChunkBuildResult& gc,
    const std::unordered_map<std::string, const GraphSite*>& site_map) {
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

        write_graph_bam_variant_one_row(out, cand, site, gc, ci);
    }
}

void write_graph_bam_variants_tsv_rows(
    std::ostream& out,
    const GraphBamChunkBuildResult& gc,
    const std::unordered_map<std::string, GraphSite>& site_by_id) {
    for (size_t ci = 0; ci < gc.chunk.candidates.size(); ++ci) {
        const CandidateVariant& cand = gc.chunk.candidates[ci];
        if (ci >= gc.site_ids.size()) continue;

        const std::string& raw_id = gc.site_ids[ci];
        auto it = site_by_id.find(raw_id);
        if (it == site_by_id.end()) {
            const auto colon_pos = raw_id.rfind(':');
            if (colon_pos != std::string::npos)
                it = site_by_id.find(raw_id.substr(0, colon_pos));
            if (it == site_by_id.end()) continue;
        }
        const GraphSite& site = it->second;
        write_graph_bam_variant_one_row(out, cand, site, gc, ci);
    }
}

void write_graph_bam_variants_tsv_rows(std::ostream& out,
                                       const GraphBamChunkBuildResult& gc) {
    if (gc.variant_emit_rows.size() != gc.chunk.candidates.size())
        throw std::runtime_error(
            "write_graph_bam_variants_tsv_rows: variant_emit_rows size mismatch");
    for (size_t ci = 0; ci < gc.chunk.candidates.size(); ++ci) {
        const GraphVariantEmitRow& row = gc.variant_emit_rows[ci];
        write_graph_bam_variant_tsv_line(out, gc.chunk.candidates[ci],
                                         row.chrom, row.pos, row.snarl_id, row.ref, row.alt);
    }
}

void write_graph_bam_variants_tsv(std::ostream& out,
                                  const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                                  const GraphSiteCatalog& catalog) {
    (void)catalog;
    write_graph_bam_variants_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_bam_variants_tsv_rows(out, gc);
}

// ── Shared BAM write helper ───────────────────────────────────────────────────

static void write_bam_record(samFile* out_sam, sam_hdr_t* hdr,
                             const std::string& name, int hap, hts_pos_t ps) {
    bam1_t* rec = bam_init1();
    if (!rec) throw std::runtime_error("bam_init1 failed");
    struct RecGuard { bam1_t* r; ~RecGuard() { bam_destroy1(r); } } rg{rec};

    if (bam_set_qname(rec, name.c_str()) < 0)
        throw std::runtime_error("bam_set_qname failed for: " + name);
    rec->core.flag = BAM_FUNMAP;
    rec->core.tid  = -1;
    rec->core.pos  = -1;
    rec->core.mtid = -1;
    rec->core.mpos = -1;
    rec->core.qual = 255;

    if (hap > 0) {
        const int32_t h = static_cast<int32_t>(hap);
        bam_aux_append(rec, "HP", 'i', sizeof(int32_t),
                       reinterpret_cast<const uint8_t*>(&h));
    }
    if (ps >= 0) {
        const int32_t p = static_cast<int32_t>(ps);
        bam_aux_append(rec, "PS", 'i', sizeof(int32_t),
                       reinterpret_cast<const uint8_t*>(&p));
    }
    if (sam_write1(out_sam, hdr, rec) < 0)
        throw std::runtime_error("failed to write BAM record for: " + name);
}

void flush_graph_phase_bam_after_merge(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::unordered_set<std::string>* next_chunk_qnames,
    std::unordered_set<std::string>& emitted_read_names) {
    std::vector<std::pair<std::string, PhaseReadOutputRow>> drained;
    drained.reserve(rows_by_read.size());
    for (auto it = rows_by_read.begin(); it != rows_by_read.end(); ) {
        if (next_chunk_qnames && next_chunk_qnames->count(it->first)) {
            ++it;
            continue;
        }
        std::string key        = std::move(it->first);
        PhaseReadOutputRow row = std::move(it->second);
        it                     = rows_by_read.erase(it);
        drained.emplace_back(std::move(key), std::move(row));
    }
    std::sort(drained.begin(), drained.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& kv : drained) {
        emitted_read_names.insert(kv.first);
        if (phase_bam_out && phase_bam_hdr) {
            const int hap =
                kv.second.has_phased_assignment ? kv.second.hap : 0;
            const hts_pos_t ps =
                kv.second.has_phased_assignment
                    ? kv.second.phase_set
                    : static_cast<hts_pos_t>(-1);
            write_bam_record(phase_bam_out, phase_bam_hdr, kv.first, hap, ps);
        }
    }
}

void write_graph_phase_bam_for_unobserved(
    samFile* phase_bam_out,
    sam_hdr_t* phase_bam_hdr,
    const std::unordered_set<std::string>& emitted_read_names,
    const std::vector<std::string>& all_read_names_sorted) {
    for (const std::string& read_name : all_read_names_sorted) {
        if (emitted_read_names.count(read_name)) continue;
        if (phase_bam_out && phase_bam_hdr)
            write_bam_record(phase_bam_out, phase_bam_hdr, read_name, 0,
                             static_cast<hts_pos_t>(-1));
    }
}

// Build sorted list of all read names: union of map keys and all_read_names.
static std::vector<std::string> sorted_read_names(
    const std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::vector<std::string>& all_read_names) {
    std::unordered_set<std::string> seen;
    seen.reserve(rows_by_read.size() + all_read_names.size());
    for (const auto& kv : rows_by_read) seen.insert(kv.first);
    for (const auto& n : all_read_names)  seen.insert(n);
    std::vector<std::string> names(seen.begin(), seen.end());
    std::sort(names.begin(), names.end());
    return names;
}

// ── End-of-run writers from pre-built rows_by_read ───────────────────────────

void write_graph_bam_phase_bam_from_rows(
    const std::string& out_path,
    const std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const std::vector<std::string>& all_read_names) {
    const auto read_names = sorted_read_names(rows_by_read, all_read_names);

    samFile* out_sam = hts_open(out_path.c_str(), "wb");
    if (!out_sam) throw std::runtime_error("failed to open output BAM: " + out_path);
    struct SamGuard { samFile* f; ~SamGuard() { hts_close(f); } } sg{out_sam};

    sam_hdr_t* hdr = sam_hdr_init();
    if (!hdr) throw std::runtime_error("failed to allocate BAM header");
    struct HdrGuard { sam_hdr_t* h; ~HdrGuard() { sam_hdr_destroy(h); } } hg{hdr};
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "queryname", nullptr);
    if (sam_hdr_write(out_sam, hdr) < 0)
        throw std::runtime_error("failed to write BAM header: " + out_path);

    for (const std::string& name : read_names) {
        auto it = rows_by_read.find(name);
        const int hap = (it != rows_by_read.end() && it->second.has_phased_assignment)
                        ? it->second.hap : 0;
        const hts_pos_t ps = (it != rows_by_read.end() && it->second.has_phased_assignment)
                             ? it->second.phase_set : static_cast<hts_pos_t>(-1);
        write_bam_record(out_sam, hdr, name, hap, ps);
    }
}

// ── Legacy multi-chunk entry points (used by collect_pipeline.cpp) ────────────

void write_graph_bam_phase_bam(const std::string& out_path,
                               const std::vector<GraphBamChunkBuildResult>& graph_chunks,
                               const std::vector<std::string>& all_read_names) {
    std::unordered_map<std::string, PhaseReadOutputRow> rows_by_read;
    for (const auto& gc : graph_chunks) merge_graph_chunk_into_read_rows(rows_by_read, gc);
    write_graph_bam_phase_bam_from_rows(out_path, rows_by_read, all_read_names);
}

} // namespace pgphase_collect
