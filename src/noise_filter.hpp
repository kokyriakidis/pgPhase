#ifndef NOISE_FILTER_HPP
#define NOISE_FILTER_HPP

/// @file noise_filter.hpp
/// @brief Reference-based noise detection for variant sites.
///
/// Shared between the BAM and graph pipelines.  All functions operate
/// on a reference sequence slice and variant coordinates — no BAM or
/// graph-specific types required.

#include <htslib/hts.h>
#include <string>
#include <vector>

#include "phasing_types.hpp"  // Interval, kSdustThreshold, kSdustWindow

namespace pgphase_collect {

/// @brief Run SDUST low-complexity detection on a reference slice.
/// @param ref_seq  Reference sequence (uppercase ACGTN).
/// @param ref_beg  1-based start of ref_seq on the chromosome.
/// @return Sorted vector of low-complexity intervals (1-based inclusive).
std::vector<Interval> find_low_complexity_intervals(
    const std::string& ref_seq,
    hts_pos_t ref_beg);

/// @brief Check if a position falls inside any low-complexity interval.
bool pos_in_low_complexity(
    hts_pos_t pos,
    const std::vector<Interval>& lc_intervals);

/// @brief Check if an indel variant sits in a homopolymer context.
///
/// Returns true if the indel's inserted/deleted bases match a
/// homopolymer run of length >= 3 at the variant position.
/// Only applied to indels with span <= @p max_xgaps.
///
/// @param pos      1-based variant position.
/// @param ref      Reference allele string.
/// @param alt      Alternate allele string.
/// @param ref_seq  Reference sequence slice.
/// @param ref_beg  1-based start of ref_seq.
/// @param ref_end  1-based end of ref_seq (inclusive).
/// @param max_xgaps  Maximum indel span to check (default 5).
bool is_homopolymer_indel(
    hts_pos_t pos,
    const std::string& ref,
    const std::string& alt,
    const std::string& ref_seq,
    hts_pos_t ref_beg,
    hts_pos_t ref_end,
    int max_xgaps = kDefaultNoisyRegMaxXgaps);

/// @brief Check if an indel variant sits in a tandem repeat context.
///
/// Returns true if the indel's inserted/deleted bases match a repeat
/// motif (1-6 bp unit) that extends >= 2 copies in the flanking
/// reference.  Only applied to indels with span <= @p max_xgaps.
///
/// @param pos      1-based variant position.
/// @param ref      Reference allele string.
/// @param alt      Alternate allele string.
/// @param ref_seq  Reference sequence slice.
/// @param ref_beg  1-based start of ref_seq.
/// @param ref_end  1-based end of ref_seq (inclusive).
/// @param max_xgaps  Maximum indel span to check (default 5).
bool is_repeat_indel(
    hts_pos_t pos,
    const std::string& ref,
    const std::string& alt,
    const std::string& ref_seq,
    hts_pos_t ref_beg,
    hts_pos_t ref_end,
    int max_xgaps = kDefaultNoisyRegMaxXgaps);

/// @brief Check if a variant site is noisy (low-complexity, homopolymer,
///        or tandem repeat context).
///
/// Combines all three checks.  For SNPs, only low-complexity is checked.
/// For indels, all three are checked.
bool is_noisy_site(
    hts_pos_t pos,
    const std::string& ref,
    const std::string& alt,
    const std::string& ref_seq,
    hts_pos_t ref_beg,
    hts_pos_t ref_end,
    const std::vector<Interval>& lc_intervals,
    int max_xgaps = kDefaultNoisyRegMaxXgaps);

} // namespace pgphase_collect

#endif // NOISE_FILTER_HPP
