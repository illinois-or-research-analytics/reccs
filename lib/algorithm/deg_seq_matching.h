#ifndef DEG_SEQ_MATCHING_H
#define DEG_SEQ_MATCHING_H

#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

void match_degree_sequence(
    Graph& g, 
    const std::shared_ptr<const std::vector<uint32_t>>& degree_sequence) {
    
    // Verify that degree_sequence is valid
    if (!degree_sequence || degree_sequence->empty()) {
        std::cerr << "Error: Degree sequence is null or empty." << std::endl;
        return;
    }

    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes in graph." 
                  << "Expected: " << g.num_nodes
                  << ", Got: " << degree_sequence->size() << std::endl;
        return;
    }

    // Random number generator for candidate selection
    std::random_device rd;
    std::mt19937 gen(rd());

    // Initialize available node degrees: {node : remaining_degree}
    std::unordered_map<uint32_t, uint32_t> available_node_degrees;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        uint32_t current_deg = g.get_degree(i);
        uint32_t target_deg = (*degree_sequence)[i];
        
        if (target_deg > current_deg) {
            available_node_degrees[i] = target_deg - current_deg;
        }
    }

    // Initialize max heap with available node degrees
    // Using negative degrees for max heap behavior with std::priority_queue
    auto heap_comparator = [](const NodeDegree& a, const NodeDegree& b) {
        if (a.degree != b.degree) return a.degree < b.degree; // Max heap
        return a.node > b.node; // Tie-breaking for deterministic behavior
    };
    
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, decltype(heap_comparator)> 
        max_heap(heap_comparator);

    // Populate initial heap
    for (const auto& pair : available_node_degrees) {
        max_heap.emplace(NodeDegree{pair.first, static_cast<int32_t>(pair.second)});
    }

    // Initialize edge lookup for fast neighbor checking
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);

    // Main algorithm loop
    while (!max_heap.empty()) {
        // Get node with highest available degree
        NodeDegree current = max_heap.top();
        max_heap.pop();
        
        uint32_t available_node = current.node;
        uint32_t avail_degree = static_cast<uint32_t>(current.degree);

        // Check if available_node is still in available_node_degrees
        if (available_node_degrees.find(available_node) == available_node_degrees.end()) {
            continue; // Node was already processed
        }

        // Verify the degree is still current (heap may contain stale entries)
        if (available_node_degrees[available_node] != avail_degree) {
            // Re-add with current degree if still positive
            if (available_node_degrees[available_node] > 0) {
                max_heap.emplace(NodeDegree{available_node, 
                    static_cast<int32_t>(available_node_degrees[available_node])});
            }
            continue;
        }

        // Find available non-neighbors
        std::vector<uint32_t> available_non_neighbors;
        available_non_neighbors.reserve(available_node_degrees.size()); // Optimize allocation
        
        for (const auto& pair : available_node_degrees) {
            uint32_t candidate = pair.first;
            if (candidate != available_node && !edge_exists_fast(available_node, candidate)) {
                available_non_neighbors.push_back(candidate);
            }
        }

        // Calculate avail_k = min(avail_degree, len(available_non_neighbors))
        uint32_t avail_k = std::min(avail_degree, static_cast<uint32_t>(available_non_neighbors.size()));

        // Connect to avail_k random non-neighbors
        std::shuffle(available_non_neighbors.begin(), available_non_neighbors.end(), gen);
        
        for (uint32_t i = 0; i < avail_k; ++i) {
            uint32_t random_node = available_non_neighbors[i];
            
            // Add edge (available_node, random_node)
            g.add_edge(available_node, random_node);
            existing_edges.insert(statics::encode_edge(available_node, random_node));

            // Update random_node's available degree
            if (available_node_degrees[random_node] > 1) {
                available_node_degrees[random_node]--;
                // Add updated node back to heap
                max_heap.emplace(NodeDegree{random_node, 
                    static_cast<int32_t>(available_node_degrees[random_node])});
            } else {
                // Delete random_node from available_node_degrees
                available_node_degrees.erase(random_node);
            }
        }

        // Delete available_node from available_node_degrees
        available_node_degrees.erase(available_node);
    }
}

#endif
