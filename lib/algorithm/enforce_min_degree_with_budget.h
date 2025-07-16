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
 * Degree-aware minimum degree enforcement
 * Replicates Python logic: prioritize available degree nodes when adding edges
 */
void enforce_min_degree_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;
    
    std::cout << "Starting simplified minimum degree enforcement on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    if (min_degree >= g.num_nodes || (g.num_nodes > 1 && min_degree > g.num_nodes - 1)) {
        std::cerr << "Error: Impossible degree requirements" << std::endl;
        return;
    }
    
    // Check if all degrees are already satisfied
    bool all_satisfied = true;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        uint32_t degree = g.row_ptr[i + 1] - g.row_ptr[i];
        if (degree < min_degree) {
            all_satisfied = false;
            break;
        }
    }
    if (all_satisfied) return;
    
    // Get cached available nodes for this cluster
    auto cluster_available_nodes = task.degree_manager->get_cluster_available_nodes(task.cluster_id);
    
    std::cout << "Cluster has " << cluster_available_nodes.size() 
              << " cached available degree nodes" << std::endl;
    
    // Pre-compute existing edges once
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track current degrees
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = g.row_ptr[i + 1] - g.row_ptr[i];
    }

    // Create priority queue for nodes needing edges
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, std::greater<NodeDegree>> degree_heap;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (current_degrees[i] < min_degree) {
            degree_heap.push({i, current_degrees[i]});
        }
    }
    
    // Pre-allocate for performance
    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
    edges_to_add.reserve(g.num_nodes * min_degree / 2);

    // Create lookup set for available nodes - one time cost
    std::unordered_set<uint64_t> available_node_set(cluster_available_nodes.begin(),
                                                   cluster_available_nodes.end());

    // Random number generator for tie-breaking
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    size_t iterations = 0;
    const size_t MAX_ITERATIONS = g.num_nodes * 3;

    // Reusable neighbor tracking
    std::vector<bool> is_neighbor(g.num_nodes, false);
    std::vector<uint32_t> neighbor_cleanup;
    neighbor_cleanup.reserve(g.num_nodes);
    
    while (!degree_heap.empty() && iterations < MAX_ITERATIONS) {
        iterations++;
        
        NodeDegree nd = degree_heap.top();
        degree_heap.pop();
        uint32_t u = nd.node;
        
        // If degree already satisfied, skip
        if (current_degrees[u] >= min_degree) continue;

        // Mark current neighbors efficiently
        neighbor_cleanup.clear();
        for (uint32_t j = g.row_ptr[u]; j < g.row_ptr[u + 1]; ++j) {
            uint32_t neighbor = g.col_idx[j];
            if (!is_neighbor[neighbor]) {
                is_neighbor[neighbor] = true;
                neighbor_cleanup.push_back(neighbor);
            }
        }
        
        // Mark self to prevent self-loops
        if (!is_neighbor[u]) {
            is_neighbor[u] = true;
            neighbor_cleanup.push_back(u);
        }
        
        // Add planned edges to neighbor set
        for (const auto& edge : edges_to_add) {
            if (edge.first == u && !is_neighbor[edge.second]) {
                is_neighbor[edge.second] = true;
                neighbor_cleanup.push_back(edge.second);
            }
            if (edge.second == u && !is_neighbor[edge.first]) {
                is_neighbor[edge.first] = true;
                neighbor_cleanup.push_back(edge.first);
            }
        }
        
        while (current_degrees[u] < min_degree) {
            // Build two lists: available non-neighbors and all non-neighbors
            std::vector<uint32_t> available_non_neighbors;
            std::vector<uint32_t> all_non_neighbors;
            
            for (uint32_t v = 0; v < g.num_nodes; ++v) {
                if (is_neighbor[v]) continue; // Skip neighbors
                
                all_non_neighbors.push_back(v);
                
                // Check if this node has available degree budget
                uint64_t v_id = g.id_map[v];
                if (available_node_set.count(v_id) && 
                    task.degree_manager->get_available_degree(v_id) > 0) {
                    available_non_neighbors.push_back(v);
                }
            }
            
            uint32_t selected_v = UINT32_MAX;
            
            // Priority 1: Choose from available non-neighbors (nodes with budget)
            if (!available_non_neighbors.empty()) {
                std::uniform_int_distribution<size_t> dist(0, available_non_neighbors.size() - 1);
                selected_v = available_non_neighbors[dist(gen)];
            }
            // Priority 2: If no available non-neighbors, choose from ALL non-neighbors
            else if (!all_non_neighbors.empty()) {
                std::uniform_int_distribution<size_t> dist(0, all_non_neighbors.size() - 1);
                selected_v = all_non_neighbors[dist(gen)];
            }
            
            if (selected_v == UINT32_MAX) {
                std::cerr << "Warning: Node " << u << " cannot reach degree " << min_degree
                          << " (no more non-neighbors available)" << std::endl;
                break; // No more candidates
            }
            
            // Add the edge
            edges_to_add.emplace_back(std::min(u, selected_v), std::max(u, selected_v));
            current_degrees[u]++;
            current_degrees[selected_v]++;
            total_edges_added++;
            
            // Track if this edge uses available degree budget
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[selected_v];
            if (available_node_set.count(u_id) || available_node_set.count(v_id)) {
                degree_corrected_edges++;
            }
            
            // Update neighbor tracking
            if (!is_neighbor[selected_v]) {
                is_neighbor[selected_v] = true;
                neighbor_cleanup.push_back(selected_v);
            }
        }
        
        // Clean up neighbor tracking
        for (uint32_t neighbor : neighbor_cleanup) {
            is_neighbor[neighbor] = false;
        }
        
        // Re-add u if still needs edges
        if (current_degrees[u] < min_degree) {
            degree_heap.push({u, current_degrees[u]});
        }
    }
    
    // Batch add all edges at once
    if (!edges_to_add.empty()) {
        add_edges_batch(g, edges_to_add);
        
        // Batch consume degree budgets efficiently
        std::vector<std::pair<uint64_t, int32_t>> consumptions;
        consumptions.reserve(edges_to_add.size() * 2);
        
        for (const auto& edge : edges_to_add) {
            uint64_t u_id = g.id_map[edge.first];
            uint64_t v_id = g.id_map[edge.second];
            
            if (available_node_set.count(u_id)) {
                consumptions.emplace_back(u_id, 1);
            }
            if (available_node_set.count(v_id)) {
                consumptions.emplace_back(v_id, 1);
            }
        }
        
        // Try batch consumption - much faster than individual attempts
        task.degree_manager->try_consume_degrees_batch(consumptions);
        
        // Update cluster cache once at the end
        task.degree_manager->update_cluster_cache(task.cluster_id, task.cluster_node_ids);
    }
    
    std::cout << "Added " << edges_to_add.size() << " edges for minimum degree " 
              << min_degree << ". Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%)" << std::endl;
}

#endif // ENFORCE_MIN_DEGREE_WITH_BUDGET_H
