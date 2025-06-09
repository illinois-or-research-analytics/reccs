#ifndef DEG_SEQ_MATCHING_H
#define DEG_SEQ_MATCHING_H

#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <set>
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

    // Initialize available node degrees map
    std::unordered_map<uint32_t, uint32_t> available_node_degrees;
    
    // Calculate remaining degrees needed for each node
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        uint32_t current_degree = g.get_degree(i);
        uint32_t target_degree = (*degree_sequence)[i];
        
        if (target_degree > current_degree) {
            available_node_degrees[i] = target_degree - current_degree;
        }
    }

    // Initialize max heap with available node degrees
    // Using negative values to simulate max heap with std::priority_queue (min heap by default)
    std::priority_queue<std::pair<int32_t, uint32_t>> max_heap;
    
    for (const auto& pair : available_node_degrees) {
        max_heap.push({static_cast<int32_t>(pair.second), pair.first}); // {-degree, node}
    }

    // Initialize edge lookup for fast checking
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Main algorithm loop
    while (!max_heap.empty()) {
        // Get node with highest available degree
        auto [avail_degree, available_node] = max_heap.top();
        max_heap.pop();
        
        // Check if node is still in available_node_degrees (may have been removed)
        if (available_node_degrees.find(available_node) == available_node_degrees.end()) {
            continue;
        }
        
        // Update avail_degree to current value (heap might be stale)
        avail_degree = available_node_degrees[available_node];
        
        // Find available non-neighbors using fast edge checking
        std::vector<uint32_t> available_non_neighbors;
        
        // Find all non-neighbors using the fast edge existence checker
        for (uint32_t node = 0; node < g.num_nodes; ++node) {
            if (node != available_node && 
                !edge_exists_fast(available_node, node) && 
                edges_to_add.count({std::min(available_node, node), std::max(available_node, node)}) == 0) {
                available_non_neighbors.push_back(node);
            }
        }
        
        // Calculate how many edges we can actually add
        uint32_t avail_k = std::min(static_cast<uint32_t>(avail_degree), 
                                   static_cast<uint32_t>(available_non_neighbors.size()));
        
        // Randomly connect to avail_k non-neighbors
        for (uint32_t i = 0; i < avail_k; ++i) {
            if (available_non_neighbors.empty()) break;
            
            // Randomly select a non-neighbor
            std::uniform_int_distribution<size_t> dist(0, available_non_neighbors.size() - 1);
            size_t random_idx = dist(gen);
            uint32_t random_node = available_non_neighbors[random_idx];
            
            // Add edge (ensure consistent ordering)
            uint32_t u = std::min(available_node, random_node);
            uint32_t v = std::max(available_node, random_node);
            edges_to_add.insert({u, v});
            
            // Remove random_node from available_non_neighbors
            available_non_neighbors.erase(available_non_neighbors.begin() + random_idx);
            
            // Update available_node_degrees for random_node
            if (available_node_degrees.find(random_node) != available_node_degrees.end()) {
                if (available_node_degrees[random_node] > 1) {
                    available_node_degrees[random_node]--;
                    // Add back to heap with updated degree
                    max_heap.push({static_cast<int32_t>(available_node_degrees[random_node]), random_node});
                } else {
                    // Delete from available_node_degrees
                    available_node_degrees.erase(random_node);
                }
            }
        }
        
        // Delete available_node from available_node_degrees
        available_node_degrees.erase(available_node);
    }

    // Add all edges in batch
    std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_vector);
    
    std::cout << "Degree sequence matching completed. Added " 
              << edges_to_add.size() << " edges." << std::endl;
}

#endif
