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

/**
 * Matches the degree sequence of a graph to a target degree sequence.
 */
void match_degree_sequence(
    Graph& g, 
    const std::shared_ptr<const std::vector<uint32_t>>& degree_sequence) {

    if (!degree_sequence || degree_sequence->empty()) {
        std::cerr << "Error: Degree sequence is null or empty." << std::endl;
        return;
    }

    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes." << std::endl;
        std::cerr << "Expected: " << g.num_nodes 
                  << ", Got: " << degree_sequence->size() << std::endl;
        return;
    }

    // Track number of edges added
    uint32_t edges_added = 0;

    // Get current degrees and sort them to match with target sequence
    std::cout << "Assembling current node degrees for graph " << g.id << "..." << std::endl;
    std::vector<std::pair<uint32_t, uint32_t>> node_degrees; // (node_id, current_degree)
    node_degrees.reserve(g.num_nodes);
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        node_degrees.emplace_back(i, g.get_degree(i));
    }
    
    // Sort by current degree in non-increasing order (same as target sequence)
    std::cout << "Sorting current node degrees for graph " << g.id << "..." << std::endl;
    std::sort(node_degrees.begin(), node_degrees.end(), 
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Compute deficits by matching sorted current degrees with target sequence
    std::cout << "Computing deficits for graph " << g.id << "..." << std::endl;
    
    // OPTIMIZATION 1: Use vectors instead of unordered_map/set for better cache performance
    std::vector<uint32_t> deficit_counts(g.num_nodes, 0);
    std::vector<uint32_t> available_nodes;
    available_nodes.reserve(g.num_nodes);
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        uint32_t node_id = node_degrees[i].first;
        uint32_t current_degree = node_degrees[i].second;
        uint32_t target_degree = (*degree_sequence)[i]; // Target sequence is already sorted
        
        if (target_degree > current_degree) {
            uint32_t deficit = target_degree - current_degree;
            deficit_counts[node_id] = deficit;
            available_nodes.push_back(node_id);
        }
    }

    if (available_nodes.empty()) {
        std::cout << "[Graph " << g.id << "]: No nodes need additional edges." << std::endl;
        return;
    }

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // OPTIMIZATION 2: Pre-allocate neighbor sets and reuse them
    std::vector<bool> is_neighbor(g.num_nodes, false);
    std::vector<uint32_t> neighbor_list; // For cleanup
    neighbor_list.reserve(1000); // Reasonable estimate
    
    std::vector<uint32_t> available_non_neighbors;
    available_non_neighbors.reserve(available_nodes.size());

    uint32_t nodes_processed = 0;

    // OPTIMIZATION 3: Process nodes in deficit order without heap
    // Sort available nodes by deficit (descending)
    std::sort(available_nodes.begin(), available_nodes.end(), 
              [&](uint32_t a, uint32_t b) { 
                  return deficit_counts[a] > deficit_counts[b]; 
              });

    // Main algorithm loop - optimized version
    std::cout << "Starting degree sequence matching for graph " << g.id << "..." << std::endl;
    
    // OPTIMIZATION 4: Use iterator-based removal instead of set operations
    auto available_end = available_nodes.end();
    
    for (auto it = available_nodes.begin(); it != available_end; ) {
        uint32_t available_node = *it;
        
        // Skip if node no longer has deficit
        if (deficit_counts[available_node] == 0) {
            ++it;
            continue;
        }
        
        uint32_t avail_degree = deficit_counts[available_node];

        // OPTIMIZATION 5: Fast neighbor marking using boolean array
        neighbor_list.clear();
        for (uint32_t idx = g.row_ptr[available_node]; idx < g.row_ptr[available_node + 1]; ++idx) {
            uint32_t neighbor = g.col_idx[idx];
            if (!is_neighbor[neighbor]) {
                is_neighbor[neighbor] = true;
                neighbor_list.push_back(neighbor);
            }
        }
        // Mark self to avoid self-loops
        if (!is_neighbor[available_node]) {
            is_neighbor[available_node] = true;
            neighbor_list.push_back(available_node);
        }

        // OPTIMIZATION 6: Build available non-neighbors efficiently
        available_non_neighbors.clear();
        for (auto it2 = available_nodes.begin(); it2 != available_end; ++it2) {
            uint32_t candidate = *it2;
            if (deficit_counts[candidate] > 0 && !is_neighbor[candidate]) {
                available_non_neighbors.push_back(candidate);
            }
        }

        uint32_t avail_k = std::min(avail_degree, static_cast<uint32_t>(available_non_neighbors.size()));
        
        // Make avail_k connections using efficient random selection
        for (uint32_t i = 0; i < avail_k && !available_non_neighbors.empty(); ++i) {
            // Use swap-and-pop for efficient random selection
            std::uniform_int_distribution<size_t> dist(0, available_non_neighbors.size() - 1);
            size_t random_index = dist(gen);
            
            uint32_t edge_end = available_non_neighbors[random_index];
            
            // Swap with last element and pop (efficient removal)
            std::swap(available_non_neighbors[random_index], available_non_neighbors.back());
            available_non_neighbors.pop_back();
            
            // Add the edge
            g.add_edge(available_node, edge_end);
            edges_added++;

            // Update edge_end's deficit
            deficit_counts[edge_end]--;
        }

        // OPTIMIZATION 7: Clean up neighbor markings efficiently
        for (uint32_t neighbor : neighbor_list) {
            is_neighbor[neighbor] = false;
        }

        // Mark this node as processed
        deficit_counts[available_node] = 0;
        
        nodes_processed++;
        if (nodes_processed % 1000 == 0) {
            // Count remaining nodes with deficits
            uint32_t remaining = 0;
            for (uint32_t node : available_nodes) {
                if (deficit_counts[node] > 0) remaining++;
            }
            std::cout << "Nodes processed: " << nodes_processed 
                      << ", Remaining nodes with deficits: " << remaining << std::endl;
        }
        
        ++it;
    }

    std::cout << "[Graph " << g.id << "]: Number of edges added: " << edges_added << std::endl;
}

#endif // DEG_SEQ_MATCHING_H
