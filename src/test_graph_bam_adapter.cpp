#include "graph_bam_adapter.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace pgphase_collect;

static bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

static GraphSite make_site(const std::string& id,
                           hts_pos_t pos,
                           const std::string& allele0,
                           const std::string& allele1) {
    GraphSite site;
    site.chrom = "chr1";
    site.ref_contig = "chr1";
    site.pos = pos;
    site.ref_beg = pos;
    site.ref_end = pos + 1;
    site.id = id;
    site.allele_traversals = {allele0, allele1};
    site.allele_walks = {parse_graph_walk(allele0), parse_graph_walk(allele1)};
    site.skip_reason = graph_site_validation_skip_reason(site);
    site.eligible = site.skip_reason.empty();
    return site;
}

int main() {
    bool ok = true;

    GraphSiteCatalog catalog;
    catalog.sites.push_back(make_site("s1", 100, ">1>2>3", ">1>4>3"));
    catalog.sites.push_back(make_site("s2", 200, ">5>6>7", ">5>8>7"));

    std::vector<GraphReadAllele> rows = {
        {"s1", "chr1", 100, "read_a", 0},
        {"s2", "chr1", 200, "read_a", 0},
        {"s1", "chr1", 100, "read_b", 0},
        {"s2", "chr1", 200, "read_b", 0},
        {"s1", "chr1", 100, "read_c", 1},
        {"s2", "chr1", 200, "read_c", 1},
        {"s1", "chr1", 100, "read_d", 1},
        {"s2", "chr1", 200, "read_d", 1},
    };

    Options build_opts;
    build_opts.min_depth = 1;
    build_opts.min_alt_depth = 1;

    std::vector<GraphChunkBuildResult> chunks;
    chunks.push_back(build_graph_chunk(catalog, rows, "chr1", 0, 300, 0, build_opts));
    ok &= check(chunks[0].chunk.candidates.size() == 2, "adapter builds two candidates");
    ok &= check(chunks[0].chunk.reads.size() == 4, "adapter builds four reads");
    ok &= check(chunks[0].chunk.read_var_profile.size() == 4, "adapter builds read profiles");
    ok &= check(chunks[0].chunk.read_var_cr != nullptr, "adapter builds read-var cgranges");

    Options opts;
    opts.read_technology = ReadTechnology::Hifi;
    phase_graph_chunks(chunks, opts);
    ok &= check(chunks[0].chunk.candidates[0].phase_set > 0, "BAM k-means phases graph candidate");
    ok &= check(chunks[0].chunk.haps.size() == 4, "BAM k-means assigns graph read hap vector");
    for (int hap : chunks[0].chunk.haps) {
        ok &= check(hap == 1 || hap == 2, "graph read hap assigned by BAM k-means");
    }

    GraphSiteCatalog conditional_catalog;
    conditional_catalog.sites.push_back(make_site("parent", 100, ">10>11>12", ">10>13>12"));
    GraphSite child = make_site("child", 110, ">20>21>22", ">20>23>22");
    child.parent = "parent";
    child.conditional_parent_alleles = {1};
    conditional_catalog.sites.push_back(child);
    std::vector<GraphReadAllele> conditional_rows = {
        {"child", "chr1", 110, "child_only", 1},
        {"parent", "chr1", 100, "with_parent", 1},
        {"child", "chr1", 110, "with_parent", 1},
    };
    Options no_af_opts;
    no_af_opts.min_depth = 1;
    no_af_opts.min_alt_depth = 1;
    no_af_opts.min_af = 0.0;
    no_af_opts.max_af = 1.0;
    GraphChunkBuildResult conditional_chunk =
        build_graph_chunk(conditional_catalog, conditional_rows, "chr1", 0, 200, 0, no_af_opts);
    ok &= check(conditional_chunk.chunk.reads.size() == 1,
                "conditional child-only graph read is missing");
    ok &= check(!conditional_chunk.chunk.reads.empty() &&
                conditional_chunk.chunk.reads[0].qname == "with_parent",
                "conditional child with parent is retained");

    // --- AF / depth filter ---
    // Four sites: one het (passes), one hom-alt (AF=1.0 > max_af), one low-depth
    // (total=4 < min_depth=5), one low-alt (alt=1 < min_alt_depth=2).
    // A cross-read observing both the het site and the hom-alt site verifies that
    // its profile is truncated to the surviving site only.
    {
        GraphSiteCatalog af_catalog;
        af_catalog.sites.push_back(make_site("af_het",       100, ">1>2>3",    ">1>4>3"));
        af_catalog.sites.push_back(make_site("af_hom_alt",   200, ">5>6>7",    ">5>8>7"));
        af_catalog.sites.push_back(make_site("af_low_depth", 300, ">9>10>11",  ">9>12>11"));
        af_catalog.sites.push_back(make_site("af_low_alt",   400, ">13>14>15", ">13>16>15"));

        std::vector<GraphReadAllele> af_rows;
        // af_het: 4 ref + 4 alt → total=8, AF=0.5 → passes
        for (int i = 0; i < 4; ++i) {
            af_rows.push_back({"af_het", "chr1", 100, "het_ref_" + std::to_string(i), 0});
            af_rows.push_back({"af_het", "chr1", 100, "het_alt_" + std::to_string(i), 1});
        }
        // af_hom_alt: 6 alt reads → AF=1.0 > max_af → filtered
        for (int i = 0; i < 6; ++i)
            af_rows.push_back({"af_hom_alt", "chr1", 200, "hom_" + std::to_string(i), 1});
        // af_low_depth: 2 ref + 2 alt → total=4 < min_depth=5 → filtered
        for (int i = 0; i < 2; ++i) {
            af_rows.push_back({"af_low_depth", "chr1", 300, "ld_ref_" + std::to_string(i), 0});
            af_rows.push_back({"af_low_depth", "chr1", 300, "ld_alt_" + std::to_string(i), 1});
        }
        // af_low_alt: 6 ref + 1 alt → alt=1 < min_alt_depth=2 → filtered
        for (int i = 0; i < 6; ++i)
            af_rows.push_back({"af_low_alt", "chr1", 400, "la_ref_" + std::to_string(i), 0});
        af_rows.push_back({"af_low_alt", "chr1", 400, "la_alt_0", 1});
        // cross_read: observes af_het (allele 0) and af_hom_alt (allele 1);
        // after filter only the af_het observation survives in its profile
        af_rows.push_back({"af_het",     "chr1", 100, "cross_read", 0});
        af_rows.push_back({"af_hom_alt", "chr1", 200, "cross_read", 1});

        Options default_opts;  // min_depth=5, min_alt_depth=2, min_af=0.20, max_af=0.80
        GraphChunkBuildResult af_chunk =
            build_graph_chunk(af_catalog, af_rows, "chr1", 0, 500, 0, default_opts);

        ok &= check(af_chunk.chunk.candidates.size() == 1,
                    "af filter: only het site survives");
        ok &= check(!af_chunk.site_ids.empty() &&
                    af_chunk.site_ids[0] == "af_het",
                    "af filter: surviving site id is af_het");
        // het_ref_0..3 + het_alt_0..3 + cross_read = 9 reads with surviving observations
        ok &= check(af_chunk.chunk.reads.size() == 9,
                    "af filter: reads with only filtered-site observations are dropped");
        ok &= check(af_chunk.chunk.candidates[0].counts.total_cov == 9,
                    "af filter: surviving site depth includes cross_read");
    }

    // --- Per-allele min_alt_depth drop ---
    // A triallelic snarl where allele 2 has only 1 supporting read (below min_alt_depth=2).
    // After dropping allele 2, the site becomes biallelic (ref=5, alt1=5) with AF=0.5 → passes.
    // Reads that observed allele 2 lose that observation and are dropped if it was their only site.
    {
        GraphSiteCatalog tri_catalog;
        // Three-allele snarl: walks ">1>2>3" (ref), ">1>4>3" (alt1), ">1>5>3" (alt2)
        GraphSite tri_site;
        tri_site.chrom = "chr1";
        tri_site.ref_contig = "chr1";
        tri_site.pos = 100;
        tri_site.ref_beg = 100;
        tri_site.ref_end = 101;
        tri_site.id = "tri";
        tri_site.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_site.allele_walks = {parse_graph_walk(">1>2>3"),
                                 parse_graph_walk(">1>4>3"),
                                 parse_graph_walk(">1>5>3")};
        tri_site.skip_reason = graph_site_validation_skip_reason(tri_site);
        tri_site.eligible = tri_site.skip_reason.empty();
        tri_catalog.sites.push_back(tri_site);

        std::vector<GraphReadAllele> tri_rows;
        // allele 0 (ref): 5 reads
        for (int i = 0; i < 5; ++i)
            tri_rows.push_back({"tri", "chr1", 100, "ref_" + std::to_string(i), 0});
        // allele 1 (alt1): 5 reads — passes min_alt_depth=2
        for (int i = 0; i < 5; ++i)
            tri_rows.push_back({"tri", "chr1", 100, "alt1_" + std::to_string(i), 1});
        // allele 2 (alt2): 1 read — below min_alt_depth=2, must be dropped
        tri_rows.push_back({"tri", "chr1", 100, "noise_read", 2});

        Options default_opts;
        GraphChunkBuildResult tri_chunk =
            build_graph_chunk(tri_catalog, tri_rows, "chr1", 0, 200, 0, default_opts);

        // Site survives: after dropping allele 2, alle_covs=[5,5], AF=0.5, total=10 >= min_depth=5
        ok &= check(tri_chunk.chunk.candidates.size() == 1,
                    "per-allele drop: triallelic site survives as biallelic after noise drop");
        ok &= check(!tri_chunk.chunk.candidates.empty() &&
                    tri_chunk.chunk.candidates[0].counts.n_uniq_alles == 2,
                    "per-allele drop: surviving site is biallelic");
        ok &= check(!tri_chunk.chunk.candidates.empty() &&
                    tri_chunk.chunk.candidates[0].counts.total_cov == 10,
                    "per-allele drop: noise read excluded from total_cov");
        // noise_read had only the dropped allele observation → no surviving site → not in reads
        ok &= check(tri_chunk.chunk.reads.size() == 10,
                    "per-allele drop: noise read with only dropped-allele observation is excluded");
    }

    // --- Multiallelic biallelic decomposition ---
    // A triallelic site (ref + alt1 + alt2) is split into two biallelic (ref vs alt_i)
    // pairs. Each pair gets its own AF/depth filter.
    // Sub-test A: both pairs pass — ref reads fan out to both pairs (allele 0 each);
    //             alt reads contribute only to their specific pair (allele 1).
    {
        GraphSite tri_both;
        tri_both.chrom = tri_both.ref_contig = "chr1";
        tri_both.pos = tri_both.ref_beg = 100; tri_both.ref_end = 101;
        tri_both.id = "tri_both";
        tri_both.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_both.allele_walks = {parse_graph_walk(">1>2>3"),
                                 parse_graph_walk(">1>4>3"),
                                 parse_graph_walk(">1>5>3")};
        tri_both.skip_reason = graph_site_validation_skip_reason(tri_both);
        tri_both.eligible = tri_both.skip_reason.empty();
        GraphSiteCatalog tri_both_cat; tri_both_cat.sites.push_back(tri_both);

        std::vector<GraphReadAllele> tb_rows;
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "ref_" + std::to_string(i), 0});
        tb_rows.push_back({"tri_both", "chr1", 100, "ref_cross", 0});
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "a1_" + std::to_string(i), 1});
        for (int i = 0; i < 5; ++i)
            tb_rows.push_back({"tri_both", "chr1", 100, "a2_" + std::to_string(i), 2});

        Options default_opts;
        auto tb = build_graph_chunk(tri_both_cat, tb_rows, "chr1", 0, 200, 0, default_opts);

        ok &= check(tb.chunk.candidates.size() == 2,
                    "multiallelic decomp: triallelic → 2 biallelic pairs");
        ok &= check(tb.site_ids.size() == 2 &&
                    tb.site_ids[0] == "tri_both:1" && tb.site_ids[1] == "tri_both:2",
                    "multiallelic decomp: pair IDs carry original alt index");
        // ref=6 (5 ref + ref_cross), alt=5 for each pair
        ok &= check(tb.chunk.candidates[0].counts.ref_cov == 6 &&
                    tb.chunk.candidates[0].counts.alt_cov == 5,
                    "multiallelic decomp: pair 0 counts correct");
        ok &= check(tb.chunk.candidates[1].counts.ref_cov == 6 &&
                    tb.chunk.candidates[1].counts.alt_cov == 5,
                    "multiallelic decomp: pair 1 counts correct");
        // All 16 reads (6 ref + 5 alt1 + 5 alt2) survive
        ok &= check(tb.chunk.reads.size() == 16,
                    "multiallelic decomp: all reads retained when both pairs pass");
    }

    // Sub-test B: alt2 has AF > max_af (hom-alt), its pair is filtered.
    // Reads observing alt2 are dropped; ref reads appear only at the surviving pair.
    {
        GraphSite tri_pf;
        tri_pf.chrom = tri_pf.ref_contig = "chr1";
        tri_pf.pos = tri_pf.ref_beg = 100; tri_pf.ref_end = 101;
        tri_pf.id = "tri_pf";
        tri_pf.allele_traversals = {">1>2>3", ">1>4>3", ">1>5>3"};
        tri_pf.allele_walks = {parse_graph_walk(">1>2>3"),
                               parse_graph_walk(">1>4>3"),
                               parse_graph_walk(">1>5>3")};
        tri_pf.skip_reason = graph_site_validation_skip_reason(tri_pf);
        tri_pf.eligible = tri_pf.skip_reason.empty();
        GraphSiteCatalog tri_pf_cat; tri_pf_cat.sites.push_back(tri_pf);

        std::vector<GraphReadAllele> pf_rows;
        for (int i = 0; i < 5; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "ref_" + std::to_string(i), 0});
        for (int i = 0; i < 5; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "a1_" + std::to_string(i), 1});
        // alt2: 50 reads → AF = 50/55 > max_af=0.80 → pair filtered
        for (int i = 0; i < 50; ++i)
            pf_rows.push_back({"tri_pf", "chr1", 100, "a2_" + std::to_string(i), 2});

        Options default_opts;
        auto pf = build_graph_chunk(tri_pf_cat, pf_rows, "chr1", 0, 200, 0, default_opts);

        ok &= check(pf.chunk.candidates.size() == 1,
                    "per-pair AF filter: only alt1 pair survives (alt2 high_af filtered)");
        ok &= check(!pf.site_ids.empty() && pf.site_ids[0] == "tri_pf:1",
                    "per-pair AF filter: surviving pair id is tri_pf:1");
        ok &= check(!pf.filtered_sites.empty() &&
                    pf.filtered_sites[0].filter_reason == "high_af",
                    "per-pair AF filter: alt2 pair recorded as high_af");
        // Only ref (5) and alt1 (5) reads survive; alt2 reads dropped
        ok &= check(pf.chunk.reads.size() == 10,
                    "per-pair AF filter: alt2 reads dropped, ref+alt1 retained");
    }

    // --- Chunk boundary: site assigned to exactly one chunk ---
    // A site whose ref_beg is at 0-based position 199 (1-based pos=200) spans the
    // boundary between chunk0=[0,200) and chunk1=[200,400).  The start-position
    // check assigns it only to chunk0.  The test pre-filters the catalog per chunk
    // (mirroring the production path) so build_graph_chunk receives only the
    // sites that belong to each chunk.
    {
        GraphSiteCatalog full_catalog;
        full_catalog.sites.push_back(make_site("left",     100, ">1>2>3",  ">1>4>3"));
        full_catalog.sites.push_back(make_site("boundary", 200, ">5>6>7",  ">5>8>7"));
        full_catalog.sites.push_back(make_site("right",    300, ">9>10>11",">9>12>11"));

        // Pre-filter catalog per chunk interval (same as production path).
        auto filter_catalog = [](const GraphSiteCatalog& cat,
                                 hts_pos_t beg, hts_pos_t end) {
            GraphSiteCatalog out;
            for (const GraphSite& site : cat.sites) {
                const hts_pos_t site_beg0 =
                    (site.ref_beg > 0 ? site.ref_beg : site.pos) - 1;
                if (site_beg0 >= beg && site_beg0 < end)
                    out.sites.push_back(site);
            }
            return out;
        };
        GraphSiteCatalog cat0 = filter_catalog(full_catalog, 0, 200);
        GraphSiteCatalog cat1 = filter_catalog(full_catalog, 200, 400);

        // Give every site 6 ref + 6 alt reads so all pass default depth/AF filters.
        std::vector<GraphReadAllele> boundary_rows;
        for (const auto& [sid, pos] : std::vector<std::pair<std::string,int>>{
                 {"left",100}, {"boundary",200}, {"right",300}}) {
            for (int i = 0; i < 6; ++i) {
                boundary_rows.push_back({sid, "chr1", pos,
                                         sid + "_ref_" + std::to_string(i), 0});
                boundary_rows.push_back({sid, "chr1", pos,
                                         sid + "_alt_" + std::to_string(i), 1});
            }
        }

        Options default_opts;
        GraphChunkBuildResult chunk0 =
            build_graph_chunk(cat0, boundary_rows, "chr1", 0,   200, 0, default_opts);
        GraphChunkBuildResult chunk1 =
            build_graph_chunk(cat1, boundary_rows, "chr1", 200, 400, 1, default_opts);

        // "left" (beg0=99) and "boundary" (beg0=199) both start in [0,200)
        ok &= check(chunk0.chunk.candidates.size() == 2,
                    "boundary: chunk0 contains left and boundary sites");
        // "right" (beg0=299) starts in [200,400); "boundary" must NOT appear again
        ok &= check(chunk1.chunk.candidates.size() == 1,
                    "boundary: chunk1 contains only right site, boundary not duplicated");
        ok &= check(!chunk1.site_ids.empty() && chunk1.site_ids[0] == "right",
                    "boundary: chunk1 site is right");
    }

    std::ostringstream sites;
    write_graph_phase_sites_tsv(sites, chunks);
    ok &= check(sites.str().find("HAP1_ALLELE") != std::string::npos,
                "BAM-derived graph phase site output");
    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
