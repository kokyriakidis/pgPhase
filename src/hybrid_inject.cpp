/// @file hybrid_inject.cpp
/// @brief Augment BAM-derived PhasingChunk with graph snarl observations.

#include <cmath>
#include "hybrid_inject.hpp"

#include "collect_phase.hpp"
#include "collect_var.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_set>

extern "C" {
#include "cgranges.h"
}

namespace pgphase_collect {

// ────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ────────────────────────────────────────────────────────────────────────────

// Declared in hybrid_inject.hpp.
VariantKey vcf_to_variant_key(int tid, hts_pos_t vcf_pos,
                              const std::string& vcf_ref,
                              const std::string& vcf_alt) {
    VariantKey key;
    key.tid = tid;
    if (vcf_ref.size() == 1 && vcf_alt.size() == 1) {
        key.type = VariantType::Snp;
        key.pos = vcf_pos;
        key.ref_len = 1;
        key.alt = vcf_alt;
    } else if (vcf_alt.size() > vcf_ref.size()) {
        // Strip the full shared prefix, mirroring the deletion branch and the
        // BAM convention (variant_key_from_digar: ref_len=0, alt=inserted
        // bases).  A clean left-anchored insertion consumes the entire REF
        // (shared == ref.size()) and yields ref_len=0; a multi-base anchor like
        // TA->TAAA reduces to pos after the shared run, alt="A", ref_len=0.  The
        // previous single-base strip mis-encoded such sites (pos and alt off by
        // the extra anchor bases).  For a single-base anchor this reduces to the
        // old result.
        size_t shared = 0;
        const size_t max_shared = std::min(vcf_ref.size(), vcf_alt.size());
        while (shared < max_shared && vcf_ref[shared] == vcf_alt[shared]) ++shared;
        key.type = VariantType::Insertion;
        key.pos = vcf_pos + static_cast<hts_pos_t>(shared);
        key.ref_len = static_cast<int>(vcf_ref.size() - shared);
        key.alt = vcf_alt.substr(shared);
    } else {
        // Strip the full shared prefix so deletions use the same normalized form
        // as the BAM path (variant_key_from_digar) and the standalone graph path
        // (graph_collect.cpp): pos = first deleted base, ref_len = deleted span,
        // alt = "" for a pure deletion.  The previous single-base anchor strip
        // left a residual base in alt and an inflated ref_len for homopolymer
        // deletions (e.g. TAA->TA yielded ref_len=2, alt="A"), which (a) blocked
        // bridging to BAM deletions in find_matching_candidate (alt never matched
        // BAM's "") and (b) collided with the BAM deletion in
        // merge_chunk_candidates, silently overwriting the verified BAM call.
        // For a single-base anchor (shared == 1) this reduces to the old result.
        size_t shared = 0;
        const size_t max_shared = std::min(vcf_ref.size(), vcf_alt.size());
        while (shared < max_shared && vcf_ref[shared] == vcf_alt[shared]) ++shared;
        key.type = VariantType::Deletion;
        key.pos = vcf_pos + static_cast<hts_pos_t>(shared);
        key.ref_len = static_cast<int>(vcf_ref.size() - shared);
        key.alt = vcf_alt.substr(shared);
    }
    return key;
}

/// Find a BAM candidate matching a VariantKey by exact position + type + allele.
///
/// The binary search requires candidates[0, n) to be sorted by sort_pos().
/// Only the original BAM candidates satisfy this (collect_var_classify sorts
/// them); graph-only candidates are appended unsorted, so callers must pass
/// n = original BAM candidate count to keep the search valid.
static int find_matching_candidate(const CandidateTable& candidates,
                                   const VariantKey& target,
                                   int n) {
    if (n > static_cast<int>(candidates.size()))
        n = static_cast<int>(candidates.size());
    const hts_pos_t target_sort = target.sort_pos();
    int lo = 0;
    int hi = n - 1;
    while (lo <= hi) {
        const int mid = lo + (hi - lo) / 2;
        const hts_pos_t mid_sort = candidates[static_cast<size_t>(mid)].key.sort_pos();
        if (mid_sort < target_sort) lo = mid + 1;
        else if (mid_sort > target_sort) hi = mid - 1;
        else {
            // Scan for exact key match at this sort_pos.
            int start = mid;
            while (start > 0 &&
                   candidates[static_cast<size_t>(start - 1)].key.sort_pos() == target_sort)
                --start;
            for (int i = start;
                 i < n &&
                 candidates[static_cast<size_t>(i)].key.sort_pos() == target_sort;
                 ++i) {
                const VariantKey& k = candidates[static_cast<size_t>(i)].key;
                if (k.type == target.type &&
                    k.ref_len == target.ref_len &&
                    k.alt == target.alt) {
                    return i;
                }
            }
            return -1;
        }
    }
    return -1;
}

/// Add a new CandidateVariant from a graph site.
///
/// The candidate is added UNCLASSIFIED (category LowCoverage, flag 0) so it is
/// excluded from k-means until its allele counts have been accumulated from
/// BAM and graph reads.  classify_graph_only_candidates() then applies the same
/// depth/AF/het gates the BAM pipeline uses before any graph site can become a
/// CleanHet phasing anchor.  Stamping CleanHet here (before counts exist) let
/// homozygous and low-support graph sites flood k-means and degrade phasing.
static int add_graph_only_candidate(PhasingChunk& chunk,
                                    const GraphSite& site,
                                    const std::string& vcf_alt,
                                    int tid) {
    const int idx = static_cast<int>(chunk.candidates.size());
    CandidateVariant cand;
    cand.key = vcf_to_variant_key(tid, site.pos, site.ref, vcf_alt);
    cand.counts.n_uniq_alles = 2;
    cand.counts.category = VariantCategory::LowCoverage;
    cand.counts.candvarcate_initial = VariantCategory::LowCoverage;
    cand.lcd_var_i_to_cate = 0;  // excluded from k-means until gated
    cand.lcd_make_variants_region_pass = true;
    chunk.candidates.push_back(std::move(cand));
    return idx;
}

// ────────────────────────────────────────────────────────────────────────────
// Phase A: site augmentation
// ────────────────────────────────────────────────────────────────────────────

SiteToCandidateMap inject_graph_sites(
        PhasingChunk& chunk,
        const GraphSiteCatalogView& graph_sites,
        const std::unordered_map<std::string, std::string>& chrom_remap,
        const Options& opts,
        int* sites_bridged_out,
        int* sites_added_out,
        std::unordered_set<int>* graph_only_candidates_out,
        GraphOnlyVcfAlleles* graph_only_vcf_alleles_out) {
    (void)opts;
    (void)chrom_remap;
    SiteToCandidateMap site_to_candidate;
    int bridged = 0, added = 0;

    // Track pre-sort indices of graph-only candidates and their original VCF
    // (ref, alt) strings so the noise filter can screen on the catalog
    // representation, matching the standalone graph pipeline.
    std::vector<int> graph_only_pre_sort;
    std::unordered_map<int, GraphOnlyVcfAllele> pre_sort_vcf_alleles;

    if (graph_sites.empty()) {
        if (sites_bridged_out) *sites_bridged_out = 0;
        if (sites_added_out) *sites_added_out = 0;
        return site_to_candidate;
    }

    const int chunk_tid = chunk.region.tid;
    const int orig_count = static_cast<int>(chunk.candidates.size());

    // find_matching_candidate binary-searches candidates[0, orig_count) and
    // requires them sorted by sort_pos().  collect_var_classify guarantees
    // this; assert it so a future change to candidate ordering fails loudly
    // here rather than silently missing bridges.
    assert(std::is_sorted(
        chunk.candidates.begin(),
        chunk.candidates.begin() + orig_count,
        [](const CandidateVariant& a, const CandidateVariant& b) {
            return a.key.sort_pos() < b.key.sort_pos();
        }));

    for (size_t si = 0; si < graph_sites.size(); ++si) {
        const GraphSite& site = graph_sites[si];
        if (!site.eligible) continue;
        if (site.alts.empty()) continue;

        // Use the same key format as GraphReadAllele.site_id.
        const std::string site_key = graph_site_key_str(site);

        // Use first non-spanning ALT allele.
        for (size_t ai = 0; ai < site.alts.size(); ++ai) {
            const std::string& vcf_alt = site.alts[ai];
            if (vcf_alt == "*") continue;

            const VariantKey target = vcf_to_variant_key(
                chunk_tid, site.pos, site.ref, vcf_alt);

            const int match_idx = find_matching_candidate(
                chunk.candidates, target, orig_count);

            if (match_idx >= 0 && match_idx < orig_count) {
                site_to_candidate[site_key] = match_idx;
                ++bridged;
            } else if (match_idx < 0) {
                const int new_idx = add_graph_only_candidate(
                    chunk, site, vcf_alt, chunk_tid);
                site_to_candidate[site_key] = new_idx;
                graph_only_pre_sort.push_back(new_idx);
                pre_sort_vcf_alleles[new_idx] =
                    GraphOnlyVcfAllele{site.pos, site.ref, vcf_alt};
                ++added;
            }
            break;  // first ALT only
        }
    }

    // Re-sort candidates if new ones were added.
    if (added > 0) {
        std::vector<size_t> order(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const hts_pos_t pa = chunk.candidates[a].key.sort_pos();
            const hts_pos_t pb = chunk.candidates[b].key.sort_pos();
            if (pa != pb) return pa < pb;
            return static_cast<int>(chunk.candidates[a].key.type) <
                   static_cast<int>(chunk.candidates[b].key.type);
        });

        std::vector<int> old_to_new(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i)
            old_to_new[order[i]] = static_cast<int>(i);

        CandidateTable sorted;
        sorted.reserve(chunk.candidates.size());
        for (size_t i = 0; i < order.size(); ++i)
            sorted.push_back(std::move(chunk.candidates[order[i]]));
        chunk.candidates = std::move(sorted);

        for (auto& [key, idx] : site_to_candidate)
            idx = old_to_new[static_cast<size_t>(idx)];

        // Remap graph-only indices through the permutation.
        if (graph_only_candidates_out) {
            for (int pre_idx : graph_only_pre_sort)
                graph_only_candidates_out->insert(
                    old_to_new[static_cast<size_t>(pre_idx)]);
        }
        if (graph_only_vcf_alleles_out) {
            for (const auto& [pre_idx, alleles] : pre_sort_vcf_alleles)
                (*graph_only_vcf_alleles_out)[
                    old_to_new[static_cast<size_t>(pre_idx)]] = alleles;
        }
    } else {
        if (graph_only_candidates_out) {
            for (int idx : graph_only_pre_sort)
                graph_only_candidates_out->insert(idx);
        }
        if (graph_only_vcf_alleles_out) {
            for (const auto& [idx, alleles] : pre_sort_vcf_alleles)
                (*graph_only_vcf_alleles_out)[idx] = alleles;
        }
    }

    if (sites_bridged_out) *sites_bridged_out = bridged;
    if (sites_added_out) *sites_added_out = added;
    return site_to_candidate;
}

// ────────────────────────────────────────────────────────────────────────────
// Phase B: graph read injection (after BAM profiles are built)
// ────────────────────────────────────────────────────────────────────────────

/// Extend an existing BAM read's profile with graph observations at
/// graph-only candidate sites.  Only fills positions where the BAM
/// profile has no informative call (allele == -1 or outside the current
/// profile span).  Returns true if any slot was filled.
static bool extend_bam_profile_with_graph_obs(
        PhasingChunk& chunk,
        int read_i,
        const std::vector<std::pair<int,int>>& graph_obs,  // (candidate_idx, allele)
        const std::unordered_set<int>& graph_only_candidates) {
    ReadVariantProfile& prof = chunk.read_var_profile[static_cast<size_t>(read_i)];

    // Track which observations are actually applied so we only update
    // allele counts for slots that were filled (not already occupied).
    std::vector<std::pair<int,int>> applied;  // (candidate_idx, allele)

    for (const auto& [cand_idx, allele] : graph_obs) {
        if (!graph_only_candidates.count(cand_idx)) continue;

        if (prof.start_var_idx < 0) {
            // Empty profile — initialize with this single observation.
            prof.start_var_idx = cand_idx;
            prof.end_var_idx = cand_idx;
            prof.alleles = {allele};
            prof.alt_qi = {0};
            applied.emplace_back(cand_idx, allele);
            continue;
        }

        if (cand_idx >= prof.start_var_idx && cand_idx <= prof.end_var_idx) {
            // Within existing span — fill only if uninformative.
            const size_t offset = static_cast<size_t>(cand_idx - prof.start_var_idx);
            if (offset < prof.alleles.size() && prof.alleles[offset] == -1) {
                prof.alleles[offset] = allele;
                applied.emplace_back(cand_idx, allele);
            }
        } else if (cand_idx < prof.start_var_idx) {
            // Extend left: prepend slots.
            const int gap = prof.start_var_idx - cand_idx;
            std::vector<int> new_alleles(static_cast<size_t>(gap), -1);
            std::vector<int> new_qi(static_cast<size_t>(gap), 0);
            new_alleles[0] = allele;
            new_alleles.insert(new_alleles.end(),
                               prof.alleles.begin(), prof.alleles.end());
            new_qi.insert(new_qi.end(),
                          prof.alt_qi.begin(), prof.alt_qi.end());
            prof.alleles = std::move(new_alleles);
            prof.alt_qi = std::move(new_qi);
            prof.start_var_idx = cand_idx;
            applied.emplace_back(cand_idx, allele);
        } else {
            // Extend right: append slots.
            const size_t new_span = static_cast<size_t>(cand_idx - prof.start_var_idx + 1);
            prof.alleles.resize(new_span, -1);
            prof.alt_qi.resize(new_span, 0);
            prof.alleles[new_span - 1] = allele;
            prof.end_var_idx = cand_idx;
            applied.emplace_back(cand_idx, allele);
        }
    }

    // Update allele counts only for observations that were actually applied.
    for (const auto& [cand_idx, allele] : applied) {
        CandidateVariant& cand =
            chunk.candidates[static_cast<size_t>(cand_idx)];
        if (allele == 0) {
            ++cand.counts.ref_cov;
            ++cand.counts.total_cov;
        } else {
            ++cand.counts.alt_cov;
            ++cand.counts.total_cov;
        }
    }

    return !applied.empty();
}

int inject_graph_reads(
        PhasingChunk& chunk,
        const std::vector<GraphReadAllele>& graph_rows,
        const SiteToCandidateMap& site_to_candidate,
        const std::unordered_set<int>& graph_only_candidates,
        const Options& opts,
        int* reads_extended_out) {
    if (graph_rows.empty() || site_to_candidate.empty()) {
        if (reads_extended_out) *reads_extended_out = 0;
        return 0;
    }

    const int chunk_tid = chunk.region.tid;

    // Build name → read index map for existing BAM reads.
    std::unordered_map<std::string, int> existing_read_idx;
    existing_read_idx.reserve(chunk.reads.size());
    for (size_t i = 0; i < chunk.reads.size(); ++i)
        existing_read_idx[chunk.reads[i].qname] = static_cast<int>(i);

    // Group graph rows by read name.
    struct ReadObs {
        int candidate_idx;
        int allele;
        int mapq;
    };
    std::unordered_map<std::string, std::vector<ReadObs>> read_observations;

    for (const GraphReadAllele& row : graph_rows) {
        if (row.allele < 0) continue;
        if (row.mapq < opts.min_mapq) continue;

        auto it = site_to_candidate.find(row.site_id);
        if (it == site_to_candidate.end()) continue;

        const int cand_idx = it->second;
        const int allele = (row.allele == 0) ? 0 : 1;
        read_observations[row.read_name].push_back(
            ReadObs{cand_idx, allele, row.mapq});
    }

    int injected = 0;
    int extended = 0;

    for (auto& [read_name, obs_vec] : read_observations) {
        if (obs_vec.empty()) continue;

        auto existing_it = existing_read_idx.find(read_name);
        if (existing_it != existing_read_idx.end()) {
            // Doubly-mapped read: extend BAM profile at graph-only sites.
            if (graph_only_candidates.empty()) continue;
            std::vector<std::pair<int,int>> graph_obs;
            graph_obs.reserve(obs_vec.size());
            for (const ReadObs& obs : obs_vec)
                graph_obs.emplace_back(obs.candidate_idx, obs.allele);
            if (extend_bam_profile_with_graph_obs(
                    chunk, existing_it->second, graph_obs,
                    graph_only_candidates)) {
                ++extended;
            }
            continue;
        }

        // Graph-only read: create synthetic read + profile.
        std::sort(obs_vec.begin(), obs_vec.end(),
                  [](const ReadObs& a, const ReadObs& b) {
                      return a.candidate_idx < b.candidate_idx;
                  });

        const int first_idx = obs_vec.front().candidate_idx;
        const int last_idx = obs_vec.back().candidate_idx;
        const int max_mapq = std::max_element(obs_vec.begin(), obs_vec.end(),
            [](const ReadObs& a, const ReadObs& b) {
                return a.mapq < b.mapq;
            })->mapq;

        const int read_id = static_cast<int>(chunk.reads.size());

        ReadRecord read;
        read.tid = chunk_tid;
        read.input_index = 0;
        read.qname = read_name;
        read.mapq = max_mapq;
        read.is_skipped = false;
        read.beg = chunk.candidates[static_cast<size_t>(first_idx)].key.sort_pos();
        read.end = chunk.candidates[static_cast<size_t>(last_idx)].key.sort_pos();
        if (read.end < read.beg) read.end = read.beg;

        ReadVariantProfile profile;
        profile.read_id = read_id;
        profile.start_var_idx = first_idx;
        profile.end_var_idx = last_idx;
        const size_t span = static_cast<size_t>(last_idx - first_idx + 1);
        profile.alleles.assign(span, -1);
        profile.alt_qi.assign(span, 0);

        for (const ReadObs& obs : obs_vec) {
            const int offset = obs.candidate_idx - first_idx;
            if (offset >= 0 && static_cast<size_t>(offset) < span)
                profile.alleles[static_cast<size_t>(offset)] = obs.allele;
        }

        chunk.reads.push_back(std::move(read));
        chunk.read_var_profile.push_back(std::move(profile));

        // Update allele counts on candidates.
        for (const ReadObs& obs : obs_vec) {
            CandidateVariant& cand =
                chunk.candidates[static_cast<size_t>(obs.candidate_idx)];
            if (obs.allele == 0) {
                ++cand.counts.ref_cov;
                ++cand.counts.total_cov;
            } else {
                ++cand.counts.alt_cov;
                ++cand.counts.total_cov;
            }
        }

        ++injected;
    }

    // Rebuild the read_var_cr interval tree with all reads (extended profiles
    // may have changed span boundaries).
    if (injected > 0 || extended > 0) {
        cgranges_t* cr = cr_init();
        if (cr == nullptr)
            throw std::runtime_error("failed to allocate hybrid read_var_cr");
        for (const ReadVariantProfile& p : chunk.read_var_profile) {
            if (p.start_var_idx < 0 || p.end_var_idx < p.start_var_idx)
                continue;
            cr_add(cr, "cr", p.start_var_idx, p.end_var_idx + 1, p.read_id);
        }
        cr_index(cr);
        chunk.read_var_cr.reset(cr);
    }

    if (reads_extended_out) *reads_extended_out = extended;
    return injected;
}

// ────────────────────────────────────────────────────────────────────────────
// Noise filter for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

void apply_hybrid_noise_filter(
        PhasingChunk& chunk,
        const std::string& ref_seq,
        hts_pos_t ref_beg,
        hts_pos_t ref_end,
        const std::unordered_set<int>& graph_only_candidates,
        int max_xgaps,
        const GraphOnlyVcfAlleles* graph_only_vcf_alleles) {
    if (ref_seq.empty() || graph_only_candidates.empty()) return;

    const std::vector<Interval> lc = find_low_complexity_intervals(ref_seq, ref_beg);

    for (int ci : graph_only_candidates) {
        if (ci < 0 || static_cast<size_t>(ci) >= chunk.candidates.size()) continue;
        CandidateVariant& cand = chunk.candidates[static_cast<size_t>(ci)];
        if (cand.counts.category != VariantCategory::CleanHetIndel) continue;

        const VariantKey& k = cand.key;
        if (k.type == VariantType::Snp) {
            // SNPs in low-complexity regions are NOT demoted: both the BAM
            // pipeline (recalls them as NOISY_CAND_HET via noisy-region MSA)
            // and the standalone graph pipeline (CLEAN_HET_SNP) keep these as
            // real het calls, so the hybrid pipeline must keep them too.
            continue;
        }

        // Prefer the original VCF (pos, ref, alt) from the snarl catalog so the
        // homopolymer/repeat verdict matches the standalone graph pipeline
        // (apply_graph_noise_filter screens on the catalog representation).
        // Reconstructing from the normalized VariantKey can yield a different
        // string/position for the same site and over-demote graph het indels
        // that the graph pipeline keeps as phasing anchors.  Fall back to
        // reconstruction only when the original is unavailable.
        hts_pos_t noisy_pos = k.pos;
        std::string vcf_ref, vcf_alt;
        bool have_alleles = false;
        if (graph_only_vcf_alleles) {
            auto it = graph_only_vcf_alleles->find(ci);
            if (it != graph_only_vcf_alleles->end()) {
                noisy_pos = it->second.pos;
                vcf_ref = it->second.ref;
                vcf_alt = it->second.alt;
                have_alleles = true;
            }
        }
        if (!have_alleles) {
            if (k.type == VariantType::Insertion) {
                const hts_pos_t anchor = k.pos - 1;
                if (anchor < ref_beg || anchor > ref_end) continue;
                const char anchor_base =
                    ref_seq[static_cast<size_t>(anchor - ref_beg)];
                vcf_ref = std::string(1, anchor_base);
                vcf_alt = anchor_base + k.alt;
            } else {
                const hts_pos_t anchor = k.pos - 1;
                if (anchor < ref_beg || k.pos + k.ref_len - 1 > ref_end) continue;
                vcf_ref = ref_seq.substr(
                    static_cast<size_t>(anchor - ref_beg),
                    static_cast<size_t>(k.ref_len + 1));
                const char anchor_base =
                    ref_seq[static_cast<size_t>(anchor - ref_beg)];
                vcf_alt = std::string(1, anchor_base);
                if (!k.alt.empty()) vcf_alt += k.alt;
            }
        }

        if (is_noisy_site(noisy_pos, vcf_ref, vcf_alt, ref_seq, ref_beg, ref_end,
                          lc, max_xgaps)) {
            cand.counts.category = VariantCategory::RepeatHetIndel;
            cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
            cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
        }
    }
}

// ────────────────────────────────────────────────────────────────────────────
// BAM count backfill for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

void backfill_graph_candidate_counts(
        PhasingChunk& chunk,
        const std::unordered_set<int>& graph_only_candidates) {
    if (graph_only_candidates.empty()) return;

    for (const ReadVariantProfile& prof : chunk.read_var_profile) {
        if (prof.start_var_idx < 0) continue;

        for (int vi = prof.start_var_idx; vi <= prof.end_var_idx; ++vi) {
            if (!graph_only_candidates.count(vi)) continue;

            const size_t offset = static_cast<size_t>(vi - prof.start_var_idx);
            if (offset >= prof.alleles.size()) continue;
            const int allele = prof.alleles[offset];

            if (allele < 0) continue;  // -1 (uninformative) or -2 (low qual)

            CandidateVariant& cand = chunk.candidates[static_cast<size_t>(vi)];
            if (allele == 0) {
                ++cand.counts.ref_cov;
            } else {
                ++cand.counts.alt_cov;
            }
            ++cand.counts.total_cov;
        }
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Quality gate for graph-only candidates
// ────────────────────────────────────────────────────────────────────────────

int classify_graph_only_candidates(
        PhasingChunk& chunk,
        const std::unordered_set<int>& graph_only_candidates,
        const Options& opts) {
    if (graph_only_candidates.empty()) return 0;

    int promoted = 0;
    for (int ci : graph_only_candidates) {
        if (ci < 0 || static_cast<size_t>(ci) >= chunk.candidates.size())
            continue;
        CandidateVariant& cand = chunk.candidates[static_cast<size_t>(ci)];
        VariantCounts& c = cand.counts;

        // Classify graph-only candidates with the SAME depth/AF/het/hom logic
        // the standalone graph pipeline uses (classify_graph_candidates), NOT
        // the BAM classify_variant_initial.  The two differ in one place: the
        // BAM classifier demotes homopolymer/repeat het indels to
        // RepeatHetIndel inline, keyed on the normalized VariantKey.  Doing
        // that here over-demotes graph het indels (the graph's main phasing
        // advantage) and pre-empts apply_hybrid_noise_filter, which screens on
        // the original catalog ref/alt.  So homopolymer screening is left
        // entirely to apply_hybrid_noise_filter (called next), mirroring the
        // graph pipeline's classify_graph_candidates + apply_graph_noise_filter
        // ordering.  Only CleanHet{Snp,Indel} enter het k-means.
        c.allele_fraction =
            c.total_cov == 0
                ? 0.0
                : static_cast<double>(c.alt_cov) /
                      static_cast<double>(c.total_cov);

        VariantCategory cat;
        if (c.total_cov < opts.min_depth || c.alt_cov < opts.min_alt_depth) {
            cat = VariantCategory::LowCoverage;
        } else if (opts.is_ont() &&
                   [&] {
                       const int fa = c.forward_alt;
                       const int ra = c.reverse_alt;
                       const int expected = (fa + ra) / 2;
                       if (expected <= 0) return false;
                       return fisher_exact_two_tail(fa, ra, expected, expected) <
                              opts.strand_bias_pval;
                   }()) {
            cat = VariantCategory::StrandBias;
        } else if (c.allele_fraction < opts.min_af) {
            // Fold LowAlleleFraction into LowCoverage so it is pruned from
            // output, matching the BAM pipeline's classify_chunk_candidates
            // pass 2 (which rewrites LOW_AF to LOW_COV before prune).
            cat = VariantCategory::LowCoverage;
        } else if (c.allele_fraction > opts.max_af) {
            cat = VariantCategory::CleanHom;
        } else if (cand.key.type == VariantType::Snp) {
            cat = VariantCategory::CleanHetSnp;
        } else if (std::abs(c.allele_fraction - 0.5) <= opts.graph_indel_af_margin &&
                   c.alt_cov >= opts.graph_indel_min_alt) {
            cat = VariantCategory::CleanHetIndel;
        } else {
            // Graph het indel with unreliable genotype (AF off-center or thin
            // alt support): keep out of k-means to avoid mis-orienting reads.
            cat = VariantCategory::LowCoverage;
        }

        c.category = cat;
        c.candvarcate_initial = cat;
        cand.lcd_var_i_to_cate = category_to_flag(cat);
        if (cat == VariantCategory::CleanHetSnp ||
            cat == VariantCategory::CleanHetIndel) {
            ++promoted;
        }
    }
    return promoted;
}

}  // namespace pgphase_collect
