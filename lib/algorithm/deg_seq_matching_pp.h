#ifndef DEG_SEQ_MATCHING_PP_H
#define DEG_SEQ_MATCHING_PP_H

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
void match_degree_sequence_pp(GraphTaskWithDegrees& task) {
    // Extract graph and degree manager from task
    Graph& g = *task.subgraph;
    auto degree_manager = task.degree_manager;
    
    if (!degree_manager) {
        std::cerr << "Error: Task's degree manager is null." << std::endl;
        return;
    }

    // Track number of edges added
    uint32_t edges_added = 0;
    uint32_t budget_consumed = 0;

    // Get available nodes for this cluster using task's helper method
    auto available_nodes = task.get_local_available_nodes();
    
    if (available_nodes.empty()) {
        std::cout << "[Cluster " << task.cluster_id << "]: No nodes need additional edges." << std::endl;
        return;
    }

    // Create local deficit map for this graph
    std::unordered_map<uint32_t, uint32_t> available_node_degrees; // local_idx -> deficit
    std::unordered_set<uint32_t> available_node_set; // local indices
    
    for (uint64_t global_node_id : available_nodes) {
        auto it = g.node_map.find(global_node_id);
        if (it != g.node_map.end()) {
            uint32_t local_idx = it->second;
            int32_t deficit = task.get_local_available_degree(global_node_id);
            if (deficit > 0) {
                available_node_degrees[local_idx] = static_cast<uint32_t>(deficit);
                available_node_set.insert(local_idx);
            }
        }
    }

    if (available_node_degrees.empty()) {
        std::cout << "[Cluster " << task.cluster_id << "]: No nodes in this cluster need additional edges." << std::endl;
        return;
    }

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // Use max heap to process largest deficits first
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

    // Main algorithm loop with budget consumption
    std::cout << "[Cluster " << task.cluster_id << "]: Starting degree deficit matching with budget tracking..." << std::endl;
    while (!max_heap.empty()) {
        NodeDegree current = max_heap.top();
        max_heap.pop();
        
        uint32_t available_node = current.node;
        uint64_t node_global_id = g.id_map[available_node];
        
        // Check if node is still available and has budget
        if (available_node_set.find(available_node) == available_node_set.end()) {
            continue;
        }
        
        // Get current available budget for this node
        int32_t current_budget = task.get_local_available_degree(node_global_id);
        if (current_budget <= 0) {
            available_node_set.erase(available_node);
            available_node_degrees.erase(available_node);
            continue;
        }
        
        uint32_t avail_degree = std::min(
            available_node_degrees[available_node], 
            static_cast<uint32_t>(current_budget)
        );

        // Get current neighbors efficiently
        std::unordered_set<uint32_t> neighbors;
        for (uint32_t idx = g.row_ptr[available_node]; idx < g.row_ptr[available_node + 1]; ++idx) {
            neighbors.insert(g.col_idx[idx]);
        }
        neighbors.insert(available_node); // Add self to avoid self-loops

        // Build available non-neighbors that also have budget
        std::vector<uint32_t> available_non_neighbors;
        available_non_neighbors.reserve(available_node_set.size());
        
        for (uint32_t candidate : available_node_set) {
            if (neighbors.find(candidate) == neighbors.end()) {
                uint64_t candidate_global_id = g.id_map[candidate];
                if (task.get_local_available_degree(candidate_global_id) > 0) {
                    available_non_neighbors.push_back(candidate);
                }
            }
        }

        uint32_t avail_k = std::min(avail_degree, static_cast<uint32_t>(available_non_neighbors.size()));
        
        // Make avail_k connections with budget consumption
        for (uint32_t i = 0; i < avail_k; ++i) {
            if (available_non_neighbors.empty()) break;
            
            // Random selection
            std::uniform_int_distribution<size_t> dist(0, available_non_neighbors.size() - 1);
            size_t idx = dist(gen);
            uint32_t edge_end = available_non_neighbors[idx];
            uint64_t edge_end_global_id = g.id_map[edge_end];
            
            // Try to consume budget for both nodes
            if (!task.consume_local_degree(node_global_id, 1)) {
                break; // No more budget for this node
            }
            
            if (!task.consume_local_degree(edge_end_global_id, 1)) {
                // Rollback the first consumption if second fails
                // Note: This is a simplified approach - in practice might need more sophisticated rollback
                continue;
            }
            
            // Successfully consumed budget, add the edge
            g.add_edge(available_node, edge_end);
            edges_added++;
            budget_consumed += 2; // Both endpoints
            
            // Swap with last and pop
            std::swap(available_non_neighbors[idx], available_non_neighbors.back());
            available_non_neighbors.pop_back();

            // Update edge_end's deficit
            available_node_degrees[edge_end]--;
            if (available_node_degrees[edge_end] == 0) {
                available_node_set.erase(edge_end);
                available_node_degrees.erase(edge_end);
            }
        }

        // Remove processed node
        available_node_set.erase(available_node);
        available_node_degrees.erase(available_node);
        
        nodes_processed++;
        if (nodes_processed % 100 == 0) {
            std::cout << "[Cluster " << task.cluster_id << "] Nodes processed: " << nodes_processed 
                      << ", Available nodes: " << available_node_set.size() << std::endl;
        }
    }

    std::cout << "[Cluster " << task.cluster_id << "]: Degree deficit matching complete. "
              << "Edges added: " << edges_added 
              << ", Budget consumed: " << budget_consumed << std::endl;
}

#endif // DEG_SEQ_MATCHING_PP_H
