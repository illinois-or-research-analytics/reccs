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

    if (!degree_sequence || degree_sequence->empty()) {
        std::cerr << "Error: Degree sequence is null or empty." << std::endl;
        return;
    }

    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes." << std::endl;
        return;
    }

    // Get current degrees and sort both sequences
    std::vector<uint32_t> current_degrees;
    current_degrees.reserve(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees.push_back(g.get_degree(i));
    }

    // Ensure both sequences are sorted – this is to minimize the number of edges added
    std::vector<uint32_t> sorted_target = *degree_sequence;
    std::vector<uint32_t> sorted_current = current_degrees;
    
    std::sort(sorted_target.begin(), sorted_target.end());
    std::sort(sorted_current.begin(), sorted_current.end());

    // Create mapping from current degrees to nodes that have those degrees
    std::unordered_map<uint32_t, std::vector<uint32_t>> degree_to_nodes;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        degree_to_nodes[g.get_degree(i)].push_back(i);
    }

    // Compute deficits using sorted sequences and assign to actual nodes
    std::unordered_map<uint32_t, uint32_t> available_node_degrees;
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        uint32_t current_sorted = sorted_current[i];
        uint32_t target_sorted = sorted_target[i];
        
        if (target_sorted > current_sorted) {
            uint32_t deficit = target_sorted - current_sorted;
            
            // Find a node with the current_sorted degree that hasn't been assigned yet
            auto& candidates = degree_to_nodes[current_sorted];
            if (!candidates.empty()) {
                uint32_t chosen_node = candidates.back();
                candidates.pop_back();
                available_node_degrees[chosen_node] = deficit;
            }
        }
    }

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // Use min heap to process smallest deficits first
    auto heap_comparator = [](const NodeDegree& a, const NodeDegree& b) {
        if (a.degree != b.degree) return a.degree > b.degree; // Min heap
        return a.node > b.node;
    };
    
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, decltype(heap_comparator)> 
        min_heap(heap_comparator);

    for (const auto& pair : available_node_degrees) {
        min_heap.emplace(NodeDegree{pair.first, static_cast<int32_t>(pair.second)});
    }

    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);

    // Main algorithm loop - corrected to match pseudocode exactly
    while (!min_heap.empty()) {
        NodeDegree current = min_heap.top();
        min_heap.pop();
        
        uint32_t available_node = current.node;
        uint32_t avail_degree = static_cast<uint32_t>(current.degree);

        // Check if available_node is still in available_node_degrees
        if (available_node_degrees.find(available_node) == available_node_degrees.end()) {
            continue;
        }

        // Check if the degree is still current (heap might have stale entries)
        if (available_node_degrees[available_node] != avail_degree) {
            if (available_node_degrees[available_node] > 0) {
                min_heap.emplace(NodeDegree{available_node, 
                    static_cast<int32_t>(available_node_degrees[available_node])});
            }
            continue;
        }

        // Find available non-neighbors (only nodes that still need more degree)
        std::vector<uint32_t> available_non_neighbors;
        for (const auto& pair : available_node_degrees) {
            uint32_t candidate = pair.first;
            if (candidate != available_node && !edge_exists_fast(available_node, candidate)) {
                available_non_neighbors.push_back(candidate);
            }
        }

        uint32_t avail_k = std::min(avail_degree, static_cast<uint32_t>(available_non_neighbors.size()));
        
        // Make avail_k connections
        for (uint32_t i = 0; i < avail_k; ++i) {
            // Randomly select from remaining available non-neighbors
            std::uniform_int_distribution<size_t> dist(0, available_non_neighbors.size() - 1);
            size_t random_index = dist(gen);
            uint32_t random_node = available_non_neighbors[random_index];
            
            // Remove the selected node from available_non_neighbors to avoid duplicate selection
            available_non_neighbors.erase(available_non_neighbors.begin() + random_index);
            
            // Add the edge
            g.add_edge(available_node, random_node);
            existing_edges.insert(statics::encode_edge(available_node, random_node));

            // Update random_node's deficit
            if (available_node_degrees[random_node] > 1) {
                available_node_degrees[random_node]--;
                // Add updated node back to heap
                min_heap.emplace(NodeDegree{random_node, 
                    static_cast<int32_t>(available_node_degrees[random_node])});
            } else {
                // random_node is satisfied, remove it
                available_node_degrees.erase(random_node);
            }
        }

        // available_node is now satisfied, remove it
        available_node_degrees.erase(available_node);
    }
}

#endif
