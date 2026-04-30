#ifndef PGPHASE_COLLECT_BAM_OUTPUT_HPP
#define PGPHASE_COLLECT_BAM_OUTPUT_HPP

#include "collect_types.hpp"

#include <memory>
#include <vector>

namespace pgphase_collect {

/** Called after closing phased/refined alignment output when `--refine-aln` is set (BAM/CRAM/SAM). */
void coordinate_sort_refined_alignment_file_or_throw(const Options& opts);

class PhasedAlignmentWriter {
public:
    PhasedAlignmentWriter(const Options& opts, const bam_hdr_t* header);
    ~PhasedAlignmentWriter();
    PhasedAlignmentWriter(const PhasedAlignmentWriter&) = delete;
    PhasedAlignmentWriter& operator=(const PhasedAlignmentWriter&) = delete;

    int write_chunks(const std::vector<BamChunk>& chunks);

private:
    Options opts_;
    htsFile* out_ = nullptr;
    std::unique_ptr<bam_hdr_t, HeaderDeleter> header_;
    std::vector<std::unique_ptr<SamFile>> in_bams_;
    std::vector<std::unique_ptr<bam_hdr_t, HeaderDeleter>> in_headers_;
    std::vector<std::unique_ptr<hts_idx_t, IndexDeleter>> in_indexes_;
};

} // namespace pgphase_collect

#endif
