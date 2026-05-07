#include "graph_phase.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <unordered_map>

namespace pgphase_collect {

namespace {

std::string site_key(const GraphSite& site, size_t index) {
    if (!site.id.empty() && site.id != ".") return site.id;
    return site.chrom + ":" + std::to_string(site.pos) + ":" + std::to_string(index);
}

std::string graph_site_root_id(const GraphSiteCatalog& catalog,
                               const std::unordered_map<std::string, size_t>& site_to_index,
                               const GraphSite& site,
                               size_t fallback_index) {
    if (!site.root.empty()) return site.root;
    if (site.parent.empty()) return site_key(site, fallback_index);
    std::string cur = site.parent;
    for (size_t depth = 0; depth < catalog.sites.size(); ++depth) {
        auto it = site_to_index.find(cur);
        if (it == site_to_index.end()) return cur;
        const GraphSite& parent = catalog.sites[it->second];
        if (parent.parent.empty()) return site_key(parent, it->second);
        cur = parent.parent;
    }
    return site_key(site, fallback_index);
}

void apply_conditional_parent_filters(GraphPhaseMatrix& matrix,
                                      const std::unordered_map<std::string, int>& site_to_phase_index) {
    for (GraphPhaseRead& read : matrix.reads) {
        std::unordered_map<int, int> allele_by_site;
        allele_by_site.reserve(read.observations.size());
        for (const GraphPhaseObservation& obs : read.observations) {
            allele_by_site.emplace(obs.site_index, obs.allele);
        }

        std::vector<GraphPhaseObservation> filtered;
        filtered.reserve(read.observations.size());
        for (const GraphPhaseObservation& obs : read.observations) {
            const GraphPhaseSite& site = matrix.sites[static_cast<size_t>(obs.site_index)];
            if (site.parent_site_id.empty() || site.conditional_parent_alleles.empty()) {
                filtered.push_back(obs);
                continue;
            }
            auto parent_it = site_to_phase_index.find(site.parent_site_id);
            if (parent_it == site_to_phase_index.end()) {
                continue;
            }
            auto parent_obs_it = allele_by_site.find(parent_it->second);
            if (parent_obs_it == allele_by_site.end()) {
                continue;
            }
            if (std::find(site.conditional_parent_alleles.begin(),
                          site.conditional_parent_alleles.end(),
                          parent_obs_it->second) != site.conditional_parent_alleles.end()) {
                filtered.push_back(obs);
            }
        }
        read.observations.swap(filtered);
    }
}

void recompute_site_allele_counts(GraphPhaseMatrix& matrix) {
    for (GraphPhaseSite& site : matrix.sites) {
        site.allele_counts.assign(static_cast<size_t>(site.n_alleles), 0);
    }
    for (const GraphPhaseRead& read : matrix.reads) {
        for (const GraphPhaseObservation& obs : read.observations) {
            if (obs.site_index < 0 || static_cast<size_t>(obs.site_index) >= matrix.sites.size()) continue;
            GraphPhaseSite& site = matrix.sites[static_cast<size_t>(obs.site_index)];
            if (obs.allele < 0 || obs.allele >= site.n_alleles) continue;
            ++site.allele_counts[static_cast<size_t>(obs.allele)];
        }
    }
}

} // namespace

GraphPhaseMatrix build_graph_phase_matrix(const GraphSiteCatalog& catalog,
                                          const std::vector<GraphReadAllele>& rows) {
    GraphPhaseMatrix matrix;
    matrix.sites.reserve(catalog.sites.size());
    std::unordered_map<std::string, int> site_to_phase_index;
    site_to_phase_index.reserve(catalog.sites.size());
    std::unordered_map<std::string, size_t> site_to_catalog_index;
    site_to_catalog_index.reserve(catalog.sites.size());

    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        const GraphSite& site = catalog.sites[i];
        const std::string key = site_key(site, i);
        if (site_to_phase_index.find(key) != site_to_phase_index.end()) {
            throw std::runtime_error("duplicate graph site id: " + key);
        }
        site_to_catalog_index.emplace(key, i);
    }

    for (size_t i = 0; i < catalog.sites.size(); ++i) {
        const GraphSite& site = catalog.sites[i];
        const std::string key = site_key(site, i);
        const int n_alleles = static_cast<int>(site.allele_walks.size());
        GraphPhaseSite out;
        out.site_id = key;
        out.chrom = site.ref_contig.empty() ? site.chrom : site.ref_contig;
        out.order_pos = site.order_pos();
        out.n_alleles = n_alleles;
        out.allele_counts.assign(static_cast<size_t>(n_alleles), 0);
        out.eligible = site.eligible && graph_site_validation_skip_reason(site).empty();
        out.level = site.level;
        out.parent_site_id = site.parent;
        out.root_site_id = graph_site_root_id(catalog, site_to_catalog_index, site, i);
        out.conditional_parent_alleles = site.conditional_parent_alleles;
        site_to_phase_index.emplace(key, static_cast<int>(matrix.sites.size()));
        matrix.sites.push_back(std::move(out));
    }

    std::unordered_map<std::string, int> read_to_index;
    read_to_index.reserve(rows.size());
    for (const GraphReadAllele& row : rows) {
        auto site_it = site_to_phase_index.find(row.site_id);
        if (site_it == site_to_phase_index.end()) continue;
        const int site_i = site_it->second;
        if (row.allele < 0 || row.allele >= matrix.sites[static_cast<size_t>(site_i)].n_alleles) continue;

        ++matrix.sites[static_cast<size_t>(site_i)].allele_counts[static_cast<size_t>(row.allele)];

        int read_i = -1;
        auto read_it = read_to_index.find(row.read_name);
        if (read_it == read_to_index.end()) {
            read_i = static_cast<int>(matrix.reads.size());
            read_to_index.emplace(row.read_name, read_i);
            GraphPhaseRead read;
            read.read_name = row.read_name;
            matrix.reads.push_back(std::move(read));
        } else {
            read_i = read_it->second;
        }
        matrix.reads[static_cast<size_t>(read_i)].observations.push_back(
            GraphPhaseObservation{site_i, row.allele});
    }

    for (GraphPhaseRead& read : matrix.reads) {
        std::sort(read.observations.begin(), read.observations.end(),
                  [](const GraphPhaseObservation& lhs, const GraphPhaseObservation& rhs) {
                      if (lhs.site_index != rhs.site_index) return lhs.site_index < rhs.site_index;
                      return lhs.allele < rhs.allele;
                  });
        std::vector<GraphPhaseObservation> dedup;
        dedup.reserve(read.observations.size());
        for (const GraphPhaseObservation& obs : read.observations) {
            if (!dedup.empty() && dedup.back().site_index == obs.site_index) {
                if (dedup.back().allele != obs.allele) dedup.back().allele = kPhaseAlleleAmbiguous;
                continue;
            }
            dedup.push_back(obs);
        }
        read.observations.clear();
        for (const GraphPhaseObservation& obs : dedup) {
            if (obs.allele >= 0) read.observations.push_back(obs);
        }
    }

    apply_conditional_parent_filters(matrix, site_to_phase_index);
    recompute_site_allele_counts(matrix);

    std::sort(matrix.reads.begin(), matrix.reads.end(),
              [](const GraphPhaseRead& lhs, const GraphPhaseRead& rhs) {
                  return lhs.read_name < rhs.read_name;
              });

    return matrix;
}

void write_graph_phase_site_counts_tsv(std::ostream& out,
                                       const GraphPhaseMatrix& matrix) {
    write_phase_site_counts_tsv(out, matrix);
}

void write_graph_phase_read_profiles_tsv(std::ostream& out,
                                         const GraphPhaseMatrix& matrix) {
    write_phase_read_profiles_tsv(out, matrix);
}

void phase_graph_matrix(GraphPhaseMatrix& matrix, int max_iters) {
    phase_matrix(matrix, max_iters);
}

void write_graph_phase_sites_tsv(std::ostream& out,
                                 const GraphPhaseMatrix& matrix) {
    write_phase_sites_tsv(out, matrix);
}

void write_graph_phase_reads_tsv(std::ostream& out,
                                 const GraphPhaseMatrix& matrix) {
    write_phase_reads_tsv(out, matrix);
}

} // namespace pgphase_collect
