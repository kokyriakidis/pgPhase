#include "graph_phase.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace pgphase_collect;

static bool check(bool cond, const std::string& msg) {
    if (!cond) std::cerr << "FAIL: " << msg << "\n";
    return cond;
}

int main() {
    GraphSiteCatalog catalog;
    {
        GraphSite s;
        s.chrom = "chr1";
        s.pos = 100;
        s.id = "s1";
        s.allele_walks = {parse_graph_walk(">1>2>3"), parse_graph_walk(">1>4>3")};
        catalog.sites.push_back(std::move(s));
    }
    {
        GraphSite s;
        s.chrom = "chr1";
        s.pos = 200;
        s.id = "s2";
        s.allele_walks = {parse_graph_walk(">5>6>7"), parse_graph_walk(">5>8>7")};
        catalog.sites.push_back(std::move(s));
    }
    {
        GraphSite s;
        s.chrom = "chr1";
        s.pos = 210;
        s.id = "s2_child";
        s.parent = "s2";
        s.level = 1;
        s.conditional_parent_alleles = {1};
        s.allele_walks = {parse_graph_walk(">9>10>11"), parse_graph_walk(">9>12>11")};
        catalog.sites.push_back(std::move(s));
    }

    std::vector<GraphReadAllele> rows = {
        {"s1", "chr1", 100, "read_a", 0, ">1>2>3"},
        {"s2", "chr1", 200, "read_a", 0, ">5>6>7"},
        {"s1", "chr1", 100, "read_b", 0, ">1>2>3"},
        {"s2", "chr1", 200, "read_b", 0, ">5>6>7"},
        {"s1", "chr1", 100, "read_c", 1, ">1>4>3"},
        {"s2", "chr1", 200, "read_c", 1, ">5>8>7"},
        {"s2_child", "chr1", 210, "read_c", 1, ">9>12>11"},
        {"s1", "chr1", 100, "read_d", 1, ">1>4>3"},
        {"s2", "chr1", 200, "read_d", 1, ">5>8>7"},
        {"s2_child", "chr1", 210, "read_d", 0, ">9>10>11"},
        {"s2_child", "chr1", 210, "read_a", 0, ">9>10>11"},
        {"s2_child", "chr1", 210, "read_e", 1, ">9>12>11"},
    };

    GraphPhaseMatrix matrix = build_graph_phase_matrix(catalog, rows);
    phase_graph_matrix(matrix);

    bool ok = true;
    ok &= check(matrix.sites.size() == 3, "three graph sites");
    ok &= check(matrix.reads.size() == 5, "five graph reads");
    ok &= check(matrix.sites[0].hap_to_cons_allele[1] != matrix.sites[0].hap_to_cons_allele[2],
                "site 1 phased");
    ok &= check(matrix.sites[1].hap_to_cons_allele[1] != matrix.sites[1].hap_to_cons_allele[2],
                "site 2 phased");
    for (const GraphPhaseRead& read : matrix.reads) {
        if (read.read_name == "read_e") {
            ok &= check(read.observations.empty(),
                        "conditional child observation dropped when parent is missing");
            ok &= check(read.hap == 0, "read with only dropped child observation is unassigned");
            continue;
        }
        ok &= check(read.hap == 1 || read.hap == 2, "read assigned hap");
        ok &= check(read.phase_set == 100, "read phase set from chunk first phased site");
        if (read.read_name == "read_a") {
            bool saw_child = false;
            for (const GraphPhaseObservation& obs : read.observations) {
                if (matrix.sites[static_cast<size_t>(obs.site_index)].site_id == "s2_child") {
                    saw_child = true;
                }
            }
            ok &= check(!saw_child, "conditional child observation dropped when parent allele mismatches");
        }
    }

    std::ostringstream sites;
    write_graph_phase_sites_tsv(sites, matrix);
    ok &= check(sites.str().find("HAP1_ALLELE") != std::string::npos,
                "phase sites header");

    std::ostringstream reads;
    write_graph_phase_reads_tsv(reads, matrix);
    ok &= check(reads.str().find("read_a") != std::string::npos,
                "phase reads row");

    GraphSiteCatalog split_catalog;
    {
        GraphSite s;
        s.chrom = "chr1";
        s.pos = 100;
        s.id = "left";
        s.allele_walks = {parse_graph_walk(">20>21>22"), parse_graph_walk(">20>23>22")};
        split_catalog.sites.push_back(std::move(s));
    }
    {
        GraphSite s;
        s.chrom = "chr1";
        s.pos = 500;
        s.id = "right";
        s.allele_walks = {parse_graph_walk(">30>31>32"), parse_graph_walk(">30>33>32")};
        split_catalog.sites.push_back(std::move(s));
    }
    std::vector<GraphReadAllele> split_rows = {
        {"left", "chr1", 100, "left_a", 0, ">20>21>22"},
        {"left", "chr1", 100, "left_b", 1, ">20>23>22"},
        {"right", "chr1", 500, "right_a", 0, ">30>31>32"},
        {"right", "chr1", 500, "right_b", 1, ">30>33>32"},
    };
    GraphPhaseMatrix split_matrix = build_graph_phase_matrix(split_catalog, split_rows);
    phase_graph_matrix(split_matrix);
    ok &= check(split_matrix.sites[0].phase_set != split_matrix.sites[1].phase_set,
                "disconnected sites keep separate phase sets");

    std::vector<PhaseChunk> chunks(2);
    chunks[0].chrom = "chr1";
    chunks[0].matrix.sites.resize(2);
    chunks[0].matrix.sites[0].site_id = "pre1";
    chunks[0].matrix.sites[0].phase_set = 100;
    chunks[0].matrix.sites[0].hap_to_cons_allele = {-1, 0, 1};
    chunks[0].matrix.sites[1].site_id = "pre2";
    chunks[0].matrix.sites[1].phase_set = 100;
    chunks[0].matrix.sites[1].hap_to_cons_allele = {-1, 0, 1};
    chunks[0].matrix.reads.push_back(ReadPhaseProfile{"shared", {}, 1, 100});

    chunks[1].chrom = "chr1";
    chunks[1].matrix.sites.resize(2);
    chunks[1].matrix.sites[0].site_id = "cur1";
    chunks[1].matrix.sites[0].phase_set = 300;
    chunks[1].matrix.sites[0].hap_to_cons_allele = {-1, 0, 1};
    chunks[1].matrix.sites[1].site_id = "cur2";
    chunks[1].matrix.sites[1].phase_set = 400;
    chunks[1].matrix.sites[1].hap_to_cons_allele = {-1, 0, 1};
    chunks[1].matrix.reads.push_back(ReadPhaseProfile{"shared", {}, 2, 300});

    stitch_phase_chunks(chunks);
    ok &= check(chunks[1].matrix.sites[0].phase_set == 100 &&
                chunks[1].matrix.sites[1].phase_set == 400,
                "stitch updates only overlapping current phase set");
    ok &= check(chunks[1].matrix.sites[0].hap_to_cons_allele[1] == 1 &&
                chunks[1].matrix.sites[1].hap_to_cons_allele[1] == 0,
                "stitch flips only overlapping current phase set");
    ok &= check(chunks[1].matrix.reads[0].hap == 1 &&
                chunks[1].matrix.reads[0].phase_set == 100,
                "stitch updates current chunk read");

    if (ok) {
        std::cout << "ALL PASS\n";
        return 0;
    }
    return 1;
}
