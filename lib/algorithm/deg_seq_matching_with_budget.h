#ifndef DEG_SEQ_MATCHING_WITH_BUDGET_H
#define DEG_SEQ_MATCHING_WITH_BUDGET_H

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
#include "../data_structures/available_node_degrees.h"

/**
 * Matches the degree sequence of a graph using pre-computed degree deficits.
 * This replaces the old degree sequence matching with node-to-node deficit consumption.
 */
void match_degree_sequence_with_budget(
    Graph& g, 
    std::shared_ptr<AvailableNodeDegreesManager> degree_manager) {

    if (!degree_manager) {
        std::cerr << "Error: Degree manager is null." << std::endl;
        return;
    }

    // Track number of edges added
    uint32_t edges_added = 0;

    // Get all nodes in this graph that have degree deficits
    std::unordered_set<uint64_t> graph_node_ids;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        graph_node_ids.insert(g.id_map[i]);
    }

    // Get available nodes for this graph
    auto available_nodes = degree_manager->get_available_nodes_for_cluster(graph_node_ids);
    
    if (available_nodes.empty()) {
        std::cout << "[Graph " << g.id << "]: No nodes need additional edges." << std::endl;
        return;
    }

    // Create local deficit map for this graph
    std::unordered_map<uint32_t, uint32_t> available_node_degrees; // local_idx -> deficit
    std::unordered_set<uint32_t> available_node_set; // local indices
    
    for (uint64_t global_node_id : available_nodes) {
        auto it = g.node_map.find(global_node_id);
        if (it != g.node_map.end()) {
            uint32_t local_idx = it->second;
            int32_t deficit = degree_manager->get_available_degree(global_node_id);
            if (deficit > 0) {
                available_node_degrees[local_idx] = static_cast<uint32_t>(deficit);
                available_node_set.insert(local_idx);
            }
        }
    }

    if (available_node_degrees.empty()) {
        std::cout << "[Graph " << g.id << "]: No nodes in this graph need additional edges." << std::endl;
        return;
    }

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // Use max heap to process largest deficits first (like Python version)
    auto heap_comparator = [](const NodeDegree& a, const NodeDegree& b) {
        if (a.degree != b.degree) return a.degree < b.degree; // Max heap
        return a.node > b.node;
    };
    
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, decltype(heap_comparator)> 
        max_heap(heap_comparator);

    for (const auto& pair : available_node_degrees) {
        max_heap.emplace(NodeDegree{pair.first, static_cast<int32_t>(pair.second)});
    }

    uint32_t nodes_processed = 0;

    // Main algorithm loop - matches Python logic exactly
    std::cout << "[Graph " << g.id << "]: Starting degree deficit matching..." << std::endl;
    while (!max_heap.empty() && !available_node_set.empty()) {
        NodeDegree current = max_heap.top();
        max_heap.pop();
        
        uint32_t available_node = current.node;
        
        // Check if node is still available
        if (available_node_set.find(available_node) == available_node_set.end()) {
            continue;
        }
        
        uint32_t avail_degree = available_node_degrees[available_node];

        // Get current neighbors efficiently
        std::unordered_set<uint32_t> neighbors;
        for (uint32_t idx = g.row_ptr[available_node]; idx < g.row_ptr[available_node + 1]; ++idx) {
            neighbors.insert(g.col_idx[idx]);
        }
        neighbors.insert(available_node); // Add self to avoid self-loops

        // Build available non-neighbors from available_node_set
        std::vector<uint32_t> available_non_neighbors;
        available_non_neighbors.reserve(available_node_set.size());
        
        for (uint32_t candidate : available_node_set) {
            if (neighbors.find(candidate) == neighbors.end()) {
                available_non_neighbors.push_back(candidate);
            }
        }

        uint32_t avail_k = std::min(avail_degree, static_cast<uint32_t>(available_non_neighbors.size()));
        
        // Make avail_k connections using efficient random selection
        for (uint32_t i = 0; i < avail_k; ++i) {
            if (available_non_neighbors.empty()) break;
            
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
            available_node_degrees[edge_end]--;
            if (available_node_degrees[edge_end] == 0) {
                available_node_set.erase(edge_end);
                available_node_degrees.erase(edge_end);
            } else {
                // Add updated node back to heap
                max_heap.emplace(NodeDegree{edge_end, 
                    static_cast<int32_t>(available_node_degrees[edge_end])});
            }
        }

        // Remove processed node
        available_node_set.erase(available_node);
        available_node_degrees.erase(available_node);
        
        nodes_processed++;
        if (nodes_processed % 100 == 0) {
            std::cout << "Nodes processed: " << nodes_processed 
                      << ", Available nodes: " << available_node_set.size() << std::endl;
        }
    }

    std::cout << "[Graph " << g.id << "]: Number of edges added for degree deficit matching: " << edges_added << std::endl;
}

#endif // DEG_SEQ_MATCHING_WITH_BUDGET_H
