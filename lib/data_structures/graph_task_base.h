#ifndef GRAPH_TASK_BASE_H
#define GRAPH_TASK_BASE_H

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include "graph.h"

// Base GraphTask struct that both workflows can inherit from or use directly
struct GraphTaskBase {
    // The subgraph to process
    std::shared_ptr<Graph> subgraph;
    
    // Cluster ID this subgraph belongs to
    std::string cluster_id;
    
    // Internal cluster index for efficiency
    uint32_t cluster_idx;
    
    // Minimum degree requirement for this cluster
    uint32_t min_degree_requirement;

    // Target degree sequence for this cluster (optional)
    std::shared_ptr<const std::vector<uint32_t>> target_degree_sequence;

    // Base constructor
    GraphTaskBase(std::shared_ptr<Graph> g, const std::string& cid, 
                  uint32_t cidx, uint32_t min_deg = 1,
                  std::shared_ptr<const std::vector<uint32_t>> deg_seq = nullptr)
        : subgraph(g), cluster_id(cid), cluster_idx(cidx), 
          min_degree_requirement(min_deg), target_degree_sequence(deg_seq) {}
};

#endif // GRAPH_TASK_BASE_H
