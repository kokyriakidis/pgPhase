/// @file noise_filter.cpp
/// @brief Reference-based noise detection — shared between BAM and graph pipelines.

#include "noise_filter.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "sdust.h"

namespace pgphase_collect {

// ────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ────────────────────────────────────────────────────────────────────────────

/// Maps ASCII base to 0–3 (ACGT) or 4 for anything else.
static int nt4_from_char(char ch) {
    switch (std::toupper(static_cast<unsigned char>(ch))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

/// Encoded base at absolute 1-based position within a reference slice.
static int ref_nt4_at(const std::string& ref_seq, hts_pos_t ref_beg, hts_pos_t abs_pos) {
    if (ref_seq.empty() || abs_pos < ref_beg) return 4;
    const hts_pos_t rel = abs_pos - ref_beg;
    if (rel < 0 || static_cast<size_t>(rel) >= ref_seq.size()) return 4;
    return nt4_from_char(ref_seq[static_cast<size_t>(rel)]);
}

/// Classify a VCF-style ref/alt pair into SNP, insertion, or deletion.
enum class VarKind { Snp, Insertion, Deletion };

static VarKind classify_variant(const std::string& ref, const std::string& alt) {
    if (ref.size() == 1 && alt.size() == 1) return VarKind::Snp;
    if (alt.size() > ref.size()) return VarKind::Insertion;
    return VarKind::Deletion;
}

// ────────────────────────────────────────────────────────────────────────────
// Public API
// ────────────────────────────────────────────────────────────────────────────

std::vector<Interval> find_low_complexity_intervals(
        const std::string& ref_seq,
        hts_pos_t ref_beg) {
    std::vector<Interval> result;
    if (ref_seq.empty()) return result;

    const int len = static_cast<int>(ref_seq.size());
    int n = 0;
    uint64_t* intervals = sdust(nullptr,
                                reinterpret_cast<const uint8_t*>(ref_seq.data()),
                                len, kSdustThreshold, kSdustWindow, &n);
    for (int i = 0; i < n; ++i) {
        const hts_pos_t rel_beg = static_cast<hts_pos_t>(intervals[i] >> 32);
        const hts_pos_t rel_end = static_cast<hts_pos_t>(static_cast<uint32_t>(intervals[i]));
        if (rel_end <= rel_beg) continue;
        result.push_back(
            Interval{ref_beg + rel_beg, ref_beg + rel_end - 1,
                     static_cast<int>(rel_end - rel_beg)});
    }
    std::free(intervals);
    return result;
}

bool pos_in_low_complexity(
        hts_pos_t pos,
        const std::vector<Interval>& lc_intervals) {
    for (const Interval& iv : lc_intervals) {
        if (pos >= iv.beg && pos <= iv.end) return true;
    }
    return false;
}

bool is_homopolymer_indel(
        hts_pos_t pos,
        const std::string& ref,
        const std::string& alt,
        const std::string& ref_seq,
        hts_pos_t ref_beg,
        hts_pos_t /*ref_end*/,
        int max_xgaps) {
    if (ref_seq.empty()) return false;

    const VarKind kind = classify_variant(ref, alt);
    hts_pos_t start_pos = 0;
    hts_pos_t end_pos = 0;

    if (kind == VarKind::Snp) {
        start_pos = pos - 1;
        end_pos = pos + 1;
    } else if (kind == VarKind::Insertion) {
        const int alt_len = static_cast<int>(alt.size()) - 1; // VCF: first base is anchor
        if (alt_len > max_xgaps) return false;
        // VCF anchor base sits at `pos`; inserted bases sit between pos and pos+1.
        // The flanking repeat context begins downstream at pos+1 and upstream at pos.
        start_pos = pos;
        end_pos = pos + 1;
    } else { // Deletion
        const int del_len = static_cast<int>(ref.size()) - 1; // VCF: first base is anchor
        if (del_len > max_xgaps) return false;
        // VCF anchor base sits at `pos`; deleted bases span pos+1..pos+del_len.
        // Downstream context resumes after the deletion; upstream context is the anchor.
        start_pos = pos;
        end_pos = pos + del_len + 1;
    }

    constexpr int max_unit_len = 6;
    constexpr int n_check_copy_num = 3;

    // Check downstream.
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i)
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, end_pos + i));

    int is_hp = 1;
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_hp = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, end_pos + i * rep_unit_len + j) != ref_bases[j]) {
                    is_hp = 0;
                    break;
                }
            }
            if (is_hp == 0) break;
        }
        if (is_hp) break;
    }
    if (is_hp) return true;

    // Check upstream.
    for (int i = 0; i < 6; ++i)
        ref_bases[i] = static_cast<uint8_t>(ref_nt4_at(ref_seq, ref_beg, start_pos - i));

    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_hp = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (ref_nt4_at(ref_seq, ref_beg, start_pos - i * rep_unit_len - j) != ref_bases[j]) {
                    is_hp = 0;
                    break;
                }
            }
            if (is_hp == 0) break;
        }
        if (is_hp) break;
    }
    return is_hp != 0;
}

bool is_repeat_indel(
        hts_pos_t pos,
        const std::string& ref,
        const std::string& alt,
        const std::string& ref_seq,
        hts_pos_t ref_beg,
        hts_pos_t ref_end,
        int max_xgaps) {
    if (ref_seq.empty()) return false;

    const VarKind kind = classify_variant(ref, alt);

    // VCF anchor base sits at `pos`; indel content begins at pos+1. The flanking
    // tandem-repeat context therefore starts one base past the anchor.
    const hts_pos_t content_pos = pos + 1;

    if (kind == VarKind::Deletion) {
        const int del_len = static_cast<int>(ref.size()) - 1; // VCF anchor base
        if (del_len <= 0 || del_len > max_xgaps) return false;
        const int len = del_len * 3;
        if (content_pos < ref_beg || content_pos + del_len + len > ref_end) return false;
        const size_t off  = static_cast<size_t>(content_pos - ref_beg);
        const size_t off2 = static_cast<size_t>(content_pos + del_len - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size() ||
            off2 + static_cast<size_t>(len) > ref_seq.size()) return false;
        return std::memcmp(ref_seq.data() + off,
                           ref_seq.data() + off2,
                           static_cast<size_t>(len)) == 0;
    }

    if (kind == VarKind::Insertion) {
        const int ins_len = static_cast<int>(alt.size()) - 1; // VCF anchor base
        if (ins_len <= 0 || ins_len > max_xgaps) return false;
        const int len = ins_len * 3;
        if (content_pos < ref_beg || content_pos + len > ref_end) return false;
        const size_t off = static_cast<size_t>(content_pos - ref_beg);
        if (off + static_cast<size_t>(len) > ref_seq.size()) return false;
        // Compare nt4-encoded so soft-masked (lowercase) reference still matches.
        const std::string ins_seq = alt.substr(1);
        for (int j = 0; j < len; ++j) {
            const int ref_base = ref_nt4_at(ref_seq, ref_beg, content_pos + j);
            const int alt_base = (j < ins_len)
                ? nt4_from_char(ins_seq[static_cast<size_t>(j)])
                : ref_nt4_at(ref_seq, ref_beg, content_pos + (j - ins_len));
            if (ref_base != alt_base) return false;
        }
        return true;
    }

    return false; // SNPs are not repeat indels.
}

bool is_noisy_site(
        hts_pos_t pos,
        const std::string& ref,
        const std::string& alt,
        const std::string& ref_seq,
        hts_pos_t ref_beg,
        hts_pos_t ref_end,
        const std::vector<Interval>& lc_intervals,
        int max_xgaps) {
    // Low-complexity applies to all variant types.
    if (pos_in_low_complexity(pos, lc_intervals)) return true;

    // Homopolymer and repeat checks apply to indels only.
    const VarKind kind = classify_variant(ref, alt);
    if (kind == VarKind::Snp) return false;

    if (is_homopolymer_indel(pos, ref, alt, ref_seq, ref_beg, ref_end, max_xgaps))
        return true;
    if (is_repeat_indel(pos, ref, alt, ref_seq, ref_beg, ref_end, max_xgaps))
        return true;

    return false;
}

} // namespace pgphase_collect
