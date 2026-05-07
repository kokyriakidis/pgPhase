#ifndef PGPHASE_GRAPH_PHASE_HPP
#define PGPHASE_GRAPH_PHASE_HPP

#include "graph_query.hpp"
#include "graph_sites.hpp"
#include "phase_engine.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace pgphase_collect {

using GraphPhaseObservation = PhaseObservation;
using GraphPhaseSite = PhaseSite;
using GraphPhaseRead = ReadPhaseProfile;
using GraphPhaseMatrix = PhaseMatrix;

GraphPhaseMatrix build_graph_phase_matrix(const GraphSiteCatalog& catalog,
                                          const std::vector<GraphReadAllele>& rows);

void write_graph_phase_site_counts_tsv(std::ostream& out,
                                       const GraphPhaseMatrix& matrix);

void write_graph_phase_read_profiles_tsv(std::ostream& out,
                                         const GraphPhaseMatrix& matrix);

void phase_graph_matrix(GraphPhaseMatrix& matrix, int max_iters = 10);

void write_graph_phase_sites_tsv(std::ostream& out,
                                 const GraphPhaseMatrix& matrix);

void write_graph_phase_reads_tsv(std::ostream& out,
                                 const GraphPhaseMatrix& matrix);

} // namespace pgphase_collect

#endif
