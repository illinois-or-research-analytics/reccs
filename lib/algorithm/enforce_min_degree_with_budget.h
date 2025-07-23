#ifndef ENFORCE_MIN_DEGREE_WITH_BUDGET_H
#define ENFORCE_MIN_DEGREE_WITH_BUDGET_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <set>
#include <random>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"
#include "../data_structures/available_node_degrees.h"

/**
 * Performance-optimized degree-aware minimum degree enforcement
 * Each task operates on its own local degree budgets - NO SHARED STATE!
 */
void enforce_min_degree_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;

    auto start_time = std::chrono::steady_clock::now();
    const auto MAX_TIME = std::chrono::seconds(30); // 30 second timeout
    
    std::cout << "Starting high-performance minimum degree enforcement on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    if (min_degree >= g.num_nodes || (g.num_nodes > 1 && min_degree > g.num_nodes - 1)) {
        std::cerr << "Error: Impossible degree requirements" << std::endl;
        return;
    }
    
    // Check if all degrees are already satisfied
    std::vector<uint32_t> current_degrees(g.num_nodes);
    std::vector<uint32_t> nodes_needing_edges;
    nodes_needing_edges.reserve(g.num_nodes);
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = g.row_ptr[i + 1] - g.row_ptr[i];
        if (current_degrees[i] < min_degree) {
            nodes_needing_edges.push_back(i);
        }
    }
    
    if (nodes_needing_edges.empty()) {
        std::cout << "All nodes already satisfy minimum degree requirement" << std::endl;
        return;
    }
    
    // OPTIMIZATION: Build fast lookup structures ONCE
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists = statics::create_edge_exists_checker(existing_edges);
    
    // OPTIMIZATION: Use the existing available nodes list from task
    task.initialize_local_degrees();
    const auto& available_node_ids = task.get_local_available_nodes(); // uint64_t global IDs
    
    // Convert to local subgraph indices and create fast lookup
    std::vector<bool> has_budget(g.num_nodes, false);
    std::vector<uint32_t> available_nodes;
    available_nodes.reserve(available_node_ids.size());
    
    for (uint64_t global_node_id : available_node_ids) {
        auto it = g.node_map.find(global_node_id);
        if (it != g.node_map.end()) {
            uint32_t local_idx = it->second;
            has_budget[local_idx] = true;
            available_nodes.push_back(local_idx);
        }
    }
    
    std::cout << "Nodes needing edges: " << nodes_needing_edges.size() 
              << ", Nodes with budget: " << available_nodes.size() << std::endl;
    
    // Track planned edges efficiently
    std::unordered_set<uint64_t> planned_edges;
    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
    edges_to_add.reserve(nodes_needing_edges.size() * min_degree);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    // OPTIMIZATION: Process nodes in batches to reduce overhead
    const size_t BATCH_SIZE = 1000;
    
    for (size_t batch_start = 0; batch_start < nodes_needing_edges.size(); batch_start += BATCH_SIZE) {
        size_t batch_end = std::min(batch_start + BATCH_SIZE, nodes_needing_edges.size());
        
        // Process batch of nodes
        for (size_t idx = batch_start; idx < batch_end; ++idx) {
            uint32_t u = nodes_needing_edges[idx];
            
            // Skip if node already satisfied
            if (current_degrees[u] >= min_degree) continue;
            
            // OPTIMIZATION: Pre-build neighbor set for this node
            std::vector<bool> is_neighbor(g.num_nodes, false);
            is_neighbor[u] = true; // Prevent self-loops
            
            // Mark existing neighbors
            for (uint32_t j = g.row_ptr[u]; j < g.row_ptr[u + 1]; ++j) {
                is_neighbor[g.col_idx[j]] = true;
            }
            
            // Mark planned neighbors
            for (const auto& edge : edges_to_add) {
                if (edge.first == u) is_neighbor[edge.second] = true;
                if (edge.second == u) is_neighbor[edge.first] = true;
            }
            
            // Add edges until degree requirement met
            while (current_degrees[u] < min_degree) {
                if (std::chrono::steady_clock::now() - start_time > MAX_TIME) {
                    std::cout << "Timeout reached, stopping min degree enforcement" << std::endl;
                    break;
                }

                uint32_t selected_v = UINT32_MAX;
                
                // STRATEGY 1: Try available nodes first (with budget)
                if (!available_nodes.empty()) {
                    std::vector<uint32_t> valid_available;
                    valid_available.reserve(available_nodes.size());
                    
                    for (uint32_t v : available_nodes) {
                        if (!is_neighbor[v] && has_budget[v]) {
                            valid_available.push_back(v);
                        }
                    }
                    
                    if (!valid_available.empty()) {
                        std::uniform_int_distribution<size_t> dist(0, valid_available.size() - 1);
                        selected_v = valid_available[dist(gen)];
                    }
                }
                
                // STRATEGY 2: If no available nodes, try any non-neighbor
                if (selected_v == UINT32_MAX) {
                    std::vector<uint32_t> all_candidates;
                    all_candidates.reserve(g.num_nodes);
                    
                    for (uint32_t v = 0; v < g.num_nodes; ++v) {
                        if (!is_neighbor[v]) {
                            all_candidates.push_back(v);
                        }
                    }
                    
                    if (!all_candidates.empty()) {
                        std::uniform_int_distribution<size_t> dist(0, all_candidates.size() - 1);
                        selected_v = all_candidates[dist(gen)];
                    }
                }
                
                if (selected_v == UINT32_MAX) {
                    std::cerr << "Warning: Node " << u << " cannot reach degree " << min_degree << std::endl;
                    break;
                }
                
                // Add edge
                uint32_t edge_u = std::min(u, selected_v);
                uint32_t edge_v = std::max(u, selected_v);
                
                edges_to_add.emplace_back(edge_u, edge_v);
                planned_edges.insert(statics::encode_edge(edge_u, edge_v));
                
                current_degrees[u]++;
                current_degrees[selected_v]++;
                total_edges_added++;
                
                // Track budget usage
                if (has_budget[u] || has_budget[selected_v]) {
                    degree_corrected_edges++;
                }
                
                // Update budget status - consume_local_degree handles list maintenance
                if (has_budget[u]) {
                    uint64_t u_id = g.id_map[u];
                    task.consume_local_degree(u_id, 1);
                    if (task.get_local_available_degree(u_id) == 0) {
                        has_budget[u] = false;
                    }
                }
                
                if (has_budget[selected_v]) {
                    uint64_t v_id = g.id_map[selected_v];
                    task.consume_local_degree(v_id, 1);
                    if (task.get_local_available_degree(v_id) == 0) {
                        has_budget[selected_v] = false;
                    }
                }
                
                // Update neighbor tracking
                is_neighbor[selected_v] = true;
            }
        }
    }
    
    // Batch add all edges
    if (!edges_to_add.empty()) {
        add_edges_batch(g, edges_to_add);
    }
    
    std::cout << "Added " << edges_to_add.size() << " edges for minimum degree " 
              << min_degree << ". Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%)" << std::endl;
}

#endif // ENFORCE_MIN_DEGREE_WITH_BUDGET_H
