#ifndef DEG_SEQ_MATCHING_H
#define DEG_SEQ_MATCHING_H

#include <vector>
#include <memory>
#include <iostream>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

void match_degree_sequence(
    Graph& g, 
    const std::shared_ptr<const std::vector<uint32_t>>& degree_sequence) {
    // Check if the degree sequence is valid
    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes in graph." << std::endl;
        return;
    }

    // Create a degree deficit vector
    // Get current degrees using NodeDegree struct
    std::vector<NodeDegree> node_degrees;
    node_degrees.reserve(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        node_degrees.emplace_back(NodeDegree{i, g.get_degree(i)});
    }

    // Sort by degree
    std::sort(node_degrees.begin(), node_degrees.end(), 
        [](const NodeDegree& a, const NodeDegree& b) {
            return a.degree < b.degree;
        });

    // Create sorted target degree sequence
    std::vector<uint32_t> sorted_target_degrees = *degree_sequence;
    std::sort(sorted_target_degrees.begin(), sorted_target_degrees.end());

    // Compute degree deficit vector using NodeDegree objects and create max heap
    std::vector<NodeDegree> degree_deficit;
    degree_deficit.reserve(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        int32_t deficit = sorted_target_degrees[i] - node_degrees[i].degree;
        degree_deficit.emplace_back(NodeDegree{node_degrees[i].node, deficit});
    }
    
    // Create max heap based on degree deficit (higher deficit = higher priority)
    std::make_heap(degree_deficit.begin(), degree_deficit.end(), 
        [](const NodeDegree& a, const NodeDegree& b) {
            return a.degree < b.degree;  // For max heap, use less-than comparator
        });

    // Initialize edge lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);

    // TODO: Implement the matching algorithm
}

#endif
