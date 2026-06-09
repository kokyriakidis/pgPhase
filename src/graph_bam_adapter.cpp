#include "graph_bam_adapter.hpp"

#include "collect_phase.hpp"
#include "fisher_exact.hpp"
#include "noise_filter.hpp"

#include <algorithm>
#include <numeric>
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

// A single read's observation at one candidate site: which site and which allele.
struct GraphProfileObservation {
    int site_index = -1;
    int allele = -1;
};

// True when a site was eligible and had ≥2 allele walks but the walk vectors
// have been cleared (CompactGraphSiteIndex took ownership of the walk data).
bool graph_site_has_released_walk_storage(const GraphSite& site) {
    if (!site.eligible || site.allele_walks.size() < 2) return false;
    return std::all_of(site.allele_walks.begin(), site.allele_walks.end(),
                       [](const GraphWalk& walk) { return walk.empty(); });
}

// Create a CandidateVariant from a GraphSite and append it to the chunk.
// Allele counts are initially zero; they are filled after parent-gated
// observation counting.
void add_graph_candidate(GraphChunkBuildResult& out,
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
    out.site_meta.push_back({
        site.ref_contig.empty() ? site.chrom : site.ref_contig,
        site.pos,
        site.ref,
        site.alts
    });
}

// Build a ReadRecord + ReadVariantProfile from a read's site observations
// and append them to the chunk.  The profile's allele vector spans from the
// first to last observed site index, with -1 for unobserved positions.
void add_read_profile(GraphChunkBuildResult& out,
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

// Build the cgranges interval index over read profiles so that
// k-means phasing can quickly find reads overlapping each candidate.
void rebuild_read_var_cr(PhasingChunk& chunk) {
    cgranges_t* cr = cr_init();
    if (cr == nullptr) throw std::runtime_error("failed to allocate graph read profile cgranges");
    for (const ReadVariantProfile& profile : chunk.read_var_profile) {
        if (profile.start_var_idx < 0 || profile.end_var_idx < profile.start_var_idx) continue;
        cr_add(cr, "cr", profile.start_var_idx, profile.end_var_idx + 1, profile.read_id);
    }
    cr_index(cr);
    chunk.read_var_cr.reset(cr);
}

// Record a read's allele at a site.  If the same site was already seen with
// a different allele (from an overlapping chunk), mark it conflicted (-1).
void merge_phase_read_observation(PhaseReadOutputRow& row,
                                  const std::string& site_id,
                                  int allele) {
    auto [it, inserted] = row.allele_by_site.emplace(site_id, allele);
    if (!inserted && it->second != allele) it->second = -1;
}

// Fold a chunk's hap/PS assignment into a read's running output row.
// Chunks are merged in order, so the last chunk's assignment wins — matching
// the BAM pipeline's PhasedAlignmentWriter where downstream ownership takes
// precedence for overlap reads.
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
void classify_graph_candidates(PhasingChunk& chunk, const Options& opts) {
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

void apply_graph_noise_filter(GraphChunkBuildResult& result,
                              const std::string& ref_seq,
                              hts_pos_t ref_beg,
                              hts_pos_t ref_end,
                              int max_xgaps) {
    if (ref_seq.empty()) return;

    const std::vector<Interval> lc = find_low_complexity_intervals(ref_seq, ref_beg);

    for (size_t ci = 0; ci < result.chunk.candidates.size(); ++ci) {
        CandidateVariant& cand = result.chunk.candidates[ci];
        // Only reclassify surviving het indels.
        if (cand.counts.category != VariantCategory::CleanHetIndel) continue;

        // Get VCF-style ref/alt from the parallel site_meta vector.
        if (ci >= result.site_meta.size()) continue;
        const GraphSiteMeta& meta = result.site_meta[ci];
        if (meta.alts.empty()) continue;

        // Use the first alt allele (biallelic pairs have exactly one).
        const std::string& vcf_ref = meta.ref;
        const std::string& vcf_alt = meta.alts[0];

        if (is_noisy_site(meta.pos, vcf_ref, vcf_alt, ref_seq, ref_beg, ref_end,
                          lc, max_xgaps)) {
            cand.counts.category = VariantCategory::RepeatHetIndel;
            cand.counts.candvarcate_initial = VariantCategory::RepeatHetIndel;
            cand.lcd_var_i_to_cate = kLongcalldRepHetVar;
        }
    }
}

// Main entry point: convert graph-space allele observations into a PhasingChunk.
//
// Pipeline within this function:
//   1. Enumerate eligible catalog sites → initial multi-allelic candidates
//   2. Wire parent→child snarl relationships for gating
//   3. Single-pass over rows: intern read names, collect per-read allele obs
//   4. Parent-gated observation counting (child obs only counted when the
//      read also observes a qualifying parent allele)
//   5. Phase 1: drop individual alt alleles below min_alt_depth
//   6. Phase 2: decompose multi-allelic sites into biallelic ref/alt pairs,
//      apply depth + AF filters per pair
//   7. Classify surviving candidates (het/hom/low-cov/strand-bias)
//   8. Phase 3: remap read observations to biallelic pair space, drop pruned
//   9. Build ReadVariantProfile + cgranges index for k-means
GraphChunkBuildResult build_graph_chunk(const GraphSiteCatalogView& catalog,
                                               const std::vector<GraphReadAllele>& rows,
                                               const std::string& /*contig*/,
                                               hts_pos_t beg,
                                               hts_pos_t end,
                                               int chunk_id,
                                               const Options& opts) {
    const int chunk_tid = 0;

    GraphChunkBuildResult out;
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
    for (size_t site_i = 0; site_i < catalog.size(); ++site_i) {
        const GraphSite& site = catalog[site_i];
        const std::string sid = graph_site_key_str(site);
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
        add_graph_candidate(out, site, sid, allele_counts.back(), chunk_tid);
    }

    // Wire parent→child snarl relationships.  A nested snarl's observations
    // are only counted when the read also traverses a qualifying allele of the
    // parent snarl (conditional_parent_alleles, from the VCF PA field).
    for (size_t site_i = 0; site_i < catalog.size(); ++site_i) {
        const GraphSite& site = catalog[site_i];
        if (site.parent.empty()) continue;
        const std::string key = graph_site_key_str(site);
        auto child_it = site_to_candidate.find(key);
        auto parent_it = site_to_candidate.find(site.parent);
        if (child_it == site_to_candidate.end() || parent_it == site_to_candidate.end()) continue;
        parent_candidate[static_cast<size_t>(child_it->second)] = parent_it->second;
    }

    // Intern read names into contiguous integer IDs to replace string-keyed maps.
    std::unordered_map<std::string, uint32_t> read_name_to_id;
    std::vector<std::string> read_id_to_name;
    read_name_to_id.reserve(rows.size());
    auto intern_read = [&](const std::string& name) -> uint32_t {
        auto [it, inserted] = read_name_to_id.emplace(name, static_cast<uint32_t>(read_id_to_name.size()));
        if (inserted) read_id_to_name.push_back(name);
        return it->second;
    };

    // Pack allele + reverse into a single int: low 16 bits = allele, bit 16 = reverse.
    // Conflict sentinel is -1 (allele < 0).
    constexpr int kRevBit = 0x10000;
    auto pack_allele_rev = [](int allele, bool rev) -> int {
        return allele | (rev ? kRevBit : 0);
    };
    auto unpack_allele = [](int packed) -> int { return packed & 0xFFFF; };
    auto unpack_reverse = [&](int packed) -> bool { return (packed & kRevBit) != 0; };

    // Single pass over rows: intern read names, track max mapq, and collect
    // per-read site allele observations. Replaces three separate passes.
    struct SiteAlleleEntry { int site_i; int packed; };
    // Temporary per-row (rid, site_i, packed) triples; sorted into per-read
    // vectors after the pass when n_reads is known.
    struct RowTriple { uint32_t rid; int site_i; int packed; int mapq; };
    std::vector<RowTriple> row_triples;
    row_triples.reserve(rows.size());
    for (const GraphReadAllele& row : rows) {
        const uint32_t rid = intern_read(row.read_name);
        auto it = site_to_candidate.find(row.site_id);
        if (it == site_to_candidate.end()) continue;
        const int site_i = it->second;
        CandidateVariant& candidate = out.chunk.candidates[static_cast<size_t>(site_i)];
        if (row.allele < 0 || row.allele >= candidate.counts.n_uniq_alles) continue;
        row_triples.push_back({rid, site_i, pack_allele_rev(row.allele, row.reverse), row.mapq});
    }
    const uint32_t n_reads = static_cast<uint32_t>(read_id_to_name.size());

    // Build max mapq and per-read allele vectors from the collected triples.
    std::vector<int> read_max_mapq(n_reads, 0);
    std::vector<std::vector<SiteAlleleEntry>> allele_by_read_site(n_reads);
    for (const RowTriple& t : row_triples) {
        if (t.mapq > read_max_mapq[t.rid]) read_max_mapq[t.rid] = t.mapq;
        allele_by_read_site[t.rid].push_back({t.site_i, t.packed});
    }
    // Sort each read's entries by site_i, then deduplicate (mark conflicts as -1).
    for (auto& entries : allele_by_read_site) {
        if (entries.empty()) continue;
        std::sort(entries.begin(), entries.end(),
                  [](const SiteAlleleEntry& a, const SiteAlleleEntry& b) {
                      return a.site_i < b.site_i;
                  });
        size_t out_i = 0;
        for (size_t j = 1; j < entries.size(); ++j) {
            if (entries[j].site_i == entries[out_i].site_i) {
                // Same site: mark conflict if alleles differ.
                if (entries[out_i].packed >= 0 &&
                    unpack_allele(entries[out_i].packed) != unpack_allele(entries[j].packed)) {
                    entries[out_i].packed = -1;
                }
            } else {
                entries[++out_i] = entries[j];
            }
        }
        entries.resize(out_i + 1);
    }

    // Binary search helper for sorted SiteAlleleEntry vectors.
    auto find_site_entry = [](const std::vector<SiteAlleleEntry>& v, int site_i)
        -> const SiteAlleleEntry* {
        auto it = std::lower_bound(v.begin(), v.end(), site_i,
                                   [](const SiteAlleleEntry& e, int s) { return e.site_i < s; });
        if (it != v.end() && it->site_i == site_i) return &*it;
        return nullptr;
    };

    // Build per-read observation lists, applying parent gating.
    std::vector<std::vector<GraphProfileObservation>> read_obs(n_reads);
    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        const auto& by_site = allele_by_read_site[rid];
        for (const auto& entry : by_site) {
            const int site_i = entry.site_i;
            const int packed = entry.packed;
            if (packed < 0) continue;
            const int allele = unpack_allele(packed);
            const bool rev = unpack_reverse(packed);
            const std::vector<int>& conditional = conditional_parent_alleles[static_cast<size_t>(site_i)];
            const int parent_i = parent_candidate[static_cast<size_t>(site_i)];
            if (!conditional.empty()) {
                if (parent_i < 0) continue;
                const SiteAlleleEntry* parent_obs = find_site_entry(by_site, parent_i);
                if (!parent_obs || parent_obs->packed < 0) continue;
                const int parent_allele = unpack_allele(parent_obs->packed);
                if (std::find(conditional.begin(), conditional.end(), parent_allele) ==
                    conditional.end()) {
                    continue;
                }
            }
            ++allele_counts[static_cast<size_t>(site_i)][static_cast<size_t>(allele)];
            auto& sc = rev ? rev_strand_counts[static_cast<size_t>(site_i)]
                           : fwd_strand_counts[static_cast<size_t>(site_i)];
            ++sc[static_cast<size_t>(allele)];
            read_obs[rid].push_back(GraphProfileObservation{site_i, allele});
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
    std::vector<GraphSiteMeta>    new_meta;
    std::vector<std::vector<int>> new_allele_counts;
    std::vector<std::vector<int>> new_fwd_strand;
    std::vector<std::vector<int>> new_rev_strand;
    std::vector<std::vector<int>> new_orig_idx;
    // Most sites are biallelic (1 surviving alt); reserve n_cands as a lower bound.
    new_cands.reserve(n_cands);
    new_ids.reserve(n_cands);
    new_meta.reserve(n_cands);
    new_allele_counts.reserve(n_cands);
    new_fwd_strand.reserve(n_cands);
    new_rev_strand.reserve(n_cands);
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
                const GraphSiteMeta& meta = out.site_meta[i];
                const size_t ref_len = meta.ref.size();
                const size_t alt_idx = static_cast<size_t>(orig_alt - 1);
                const size_t alt_len = alt_idx < meta.alts.size() ? meta.alts[alt_idx].size() : ref_len;
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
            new_meta.push_back(out.site_meta[i]);
            new_allele_counts.push_back({ref_c, alt_c});
            new_fwd_strand.push_back({fwd_ref, fwd_alt});
            new_rev_strand.push_back({rev_ref, rev_alt});
            new_orig_idx.push_back({0, orig_alt});
        }
    }

    out.chunk.candidates       = std::move(new_cands);
    out.site_ids               = std::move(new_ids);
    out.site_meta              = std::move(new_meta);
    allele_counts              = std::move(new_allele_counts);
    fwd_strand_counts          = std::move(new_fwd_strand);
    rev_strand_counts          = std::move(new_rev_strand);
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
    for (uint32_t rid = 0; rid < n_reads; ++rid) {
        auto& obs = read_obs[rid];
        std::vector<GraphProfileObservation> new_obs;
        new_obs.reserve(obs.size());
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

    // Build sorted read order by name for deterministic output.
    std::vector<uint32_t> read_order(n_reads);
    std::iota(read_order.begin(), read_order.end(), 0);
    std::sort(read_order.begin(), read_order.end(),
              [&](uint32_t a, uint32_t b) { return read_id_to_name[a] < read_id_to_name[b]; });

    for (uint32_t rid : read_order) {
        auto& observations = read_obs[rid];
        if (observations.empty()) continue;
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
        const int mapq = read_max_mapq[rid];
        add_read_profile(out, read_id_to_name[rid], dedup, mapq);
    }

    out.chunk.haps.assign(out.chunk.reads.size(), 0);
    out.chunk.phase_sets.assign(out.chunk.reads.size(), -1);
    rebuild_read_var_cr(out.chunk);
    out.site_allele_orig_idx = std::move(allele_orig_idx);

    return out;
}

// Merge-intersect overlap detection between adjacent chunks.
// Reads are inserted in sorted name order by build_graph_chunk,
// so we can merge-intersect directly in O(n+m) without sorting.
static void populate_graph_chunk_pair_overlap_impl(PhasingChunk& pre, PhasingChunk& cur) {
    pre.down_ovlp_read_i.assign(1, {});
    cur.up_ovlp_read_i.assign(1, {});
    pre.n_down_ovlp_skip_reads.assign(1, 0);
    cur.n_up_ovlp_skip_reads.assign(1, 0);

    size_t pi = 0, ci = 0;
    while (pi < pre.reads.size() && ci < cur.reads.size()) {
        const int cmp = pre.reads[pi].qname.compare(cur.reads[ci].qname);
        if (cmp < 0) { ++pi; }
        else if (cmp > 0) { ++ci; }
        else {
            pre.down_ovlp_read_i[0].push_back(static_cast<int>(pi));
            cur.up_ovlp_read_i[0].push_back(static_cast<int>(ci));
            ++pi; ++ci;
        }
    }
}

void populate_graph_chunk_overlaps(std::vector<GraphChunkBuildResult>& graph_chunks) {
    for (size_t i = 1; i < graph_chunks.size(); ++i)
        populate_graph_chunk_pair_overlap_impl(graph_chunks[i - 1].chunk, graph_chunks[i].chunk);
}

void phase_graph_chunks(std::vector<GraphChunkBuildResult>& graph_chunks,
                            const Options& opts) {
    for (GraphChunkBuildResult& graph_chunk : graph_chunks) {
        assign_hap_based_on_germline_het_vars_kmeans(graph_chunk.chunk, opts, kCandGermlineClean);
    }
    populate_graph_chunk_overlaps(graph_chunks);
    std::vector<PhasingChunk> chunks;
    chunks.reserve(graph_chunks.size());
    for (GraphChunkBuildResult& graph_chunk : graph_chunks) {
        chunks.push_back(std::move(graph_chunk.chunk));
    }
    stitch_chunk_haps(chunks, &opts, nullptr);
    for (size_t i = 0; i < chunks.size(); ++i) {
        graph_chunks[i].chunk = std::move(chunks[i]);
    }
}

void merge_graph_chunk_into_read_rows(
    std::unordered_map<std::string, PhaseReadOutputRow>& rows_by_read,
    const GraphChunkBuildResult& gc) {
    const PhasingChunk& chunk = gc.chunk;
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

void write_graph_phase_sites_tsv_header(std::ostream& out) {
    out << "CHUNK_ID\tSITE_INDEX\tSITE_ID\tPOS\tN_ALLELES\tDEPTH\tALLELE_COUNTS\tPHASE_SET\tHAP1_ALLELE\tHAP2_ALLELE\n";
}

void write_graph_phase_sites_tsv_rows(std::ostream& out,
                                          const GraphChunkBuildResult& gc) {
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

void write_graph_phase_sites_tsv(std::ostream& out,
                                     const std::vector<GraphChunkBuildResult>& graph_chunks) {
    write_graph_phase_sites_tsv_header(out);
    for (const auto& gc : graph_chunks)
        write_graph_phase_sites_tsv_rows(out, gc);
}

// Write a minimal unmapped BAM record carrying HP and PS aux tags.
// Used by flush_graph_phase_bam_after_merge to emit phased read assignments.
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

} // namespace pgphase_collect
