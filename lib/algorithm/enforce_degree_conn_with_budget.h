#ifndef ENFORCE_DEGREE_CONN_WITH_BUDGET_H
#define ENFORCE_DEGREE_CONN_WITH_BUDGET_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <random>
#include <set>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"
#include "../data_structures/available_node_degrees.h"
#include "../utils/statics.h"
#include "pcg_random.hpp"

// Helper: Find connected components using BFS
std::vector<std::vector<uint32_t>> find_connected_components_budget(const Graph& g) {
    std::vector<bool> visited(g.num_nodes, false);
    std::vector<std::vector<uint32_t>> components;
    
    for (uint32_t start = 0; start < g.num_nodes; ++start) {
        if (visited[start]) continue;
        
        std::vector<uint32_t> component;
        std::queue<uint32_t> q;
        q.push(start);
        visited[start] = true;
        
        while (!q.empty()) {
            uint32_t u = q.front();
            q.pop();
            component.push_back(u);
            
            for (uint32_t i = g.row_ptr[u]; i < g.row_ptr[u + 1]; ++i) {
                uint32_t v = g.col_idx[i];
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        
        components.push_back(std::move(component));
    }
    
    return components;
}

/**
 * Combined degree enforcement and connectivity with budget awareness
 * Uses randomized edge selection (like reference files) to avoid inflating diameter
 */
void enforce_degree_and_connectivity_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;
    
    auto start_time = std::chrono::steady_clock::now();
    const auto MAX_TIME = std::chrono::seconds(3600);
    
    std::cout << "[Cluster " << task.cluster_id << "] Starting budget-aware degree and connectivity enforcement "
              << "(min_degree=" << min_degree << ")" << std::endl;
    
    // Check if the requirement is even possible
    if (min_degree >= g.num_nodes) {
        std::cerr << "[Cluster " << task.cluster_id << "] Error: min_degree " << min_degree 
                  << " >= num_nodes " << g.num_nodes << std::endl;
        return;
    }
    
    if (g.num_nodes > 1 && min_degree > g.num_nodes - 1) {
        std::cerr << "[Cluster " << task.cluster_id << "] Error: Impossible to achieve min_degree " << min_degree 
                  << " with " << g.num_nodes << " nodes" << std::endl;
        return;
    }
    
    // Get nodes with available budget
    const auto& available_node_ids = task.get_local_available_nodes();
    
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
    
    std::cout << "[Cluster " << task.cluster_id << "] " << available_nodes.size() 
              << " nodes have degree budget available" << std::endl;
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add
    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
    std::unordered_set<uint64_t> planned_edges;
    
    // Track current degrees (will update as we add edges)
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = g.get_degree(i);
    }
    
    // Random number generator
    std::random_device rd;
    pcg32 gen(rd());
    
    // Track statistics
    uint32_t connectivity_edges = 0;
    uint32_t degree_edges = 0;
    uint32_t budget_consumed = 0;
    
    // ==================================================================
    // PHASE 1: CONNECTIVITY ENFORCEMENT (Keep same as original)
    // ==================================================================
    
    bool timeout_reached = false;
    
    auto components = find_connected_components_budget(g);
    
    std::cout << "[Cluster " << task.cluster_id << "] Found " << components.size() 
              << " connected components" << std::endl;
    
    if (components.size() > 1) {
        // Sort components by size
        std::sort(components.begin(), components.end(), 
                  [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                      return a.size() > b.size();
                  });
        
        // Connect all components to the largest component
        for (size_t i = 1; i < components.size() && !timeout_reached; ++i) {
            const auto& component_a = components[0];
            const auto& component_b = components[i];
            
            // Draw min_degree edges between components
            for (uint32_t edge_count = 0; edge_count < min_degree; ++edge_count) {
                if (std::chrono::steady_clock::now() - start_time > MAX_TIME) {
                    std::cout << "[Cluster " << task.cluster_id << "] Timeout during connectivity" << std::endl;
                    timeout_reached = true;
                    break;
                }
                
                // Build candidate lists
                std::vector<uint32_t> available_nodes_a, all_nodes_a;
                std::vector<uint32_t> available_nodes_b, all_nodes_b;
                
                for (uint32_t node : component_a) {
                    all_nodes_a.push_back(node);
                    if (has_budget[node]) {
                        available_nodes_a.push_back(node);
                    }
                }
                
                for (uint32_t node : component_b) {
                    all_nodes_b.push_back(node);
                    if (has_budget[node]) {
                        available_nodes_b.push_back(node);
                    }
                }
                
                // Select first endpoint
                uint32_t selected_a;
                bool a_has_budget = false;
                
                if (!available_nodes_a.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, available_nodes_a.size() - 1);
                    selected_a = available_nodes_a[dist(gen)];
                    a_has_budget = true;
                } else if (!all_nodes_a.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, all_nodes_a.size() - 1);
                    selected_a = all_nodes_a[dist(gen)];
                } else {
                    break;
                }
                
                // Filter component B to exclude neighbors
                std::unordered_set<uint32_t> neighbors_of_a;
                for (uint32_t idx = g.row_ptr[selected_a]; idx < g.row_ptr[selected_a + 1]; ++idx) {
                    neighbors_of_a.insert(g.col_idx[idx]);
                }
                
                std::vector<uint32_t> filtered_available_b, filtered_all_b;
                for (uint32_t node : available_nodes_b) {
                    if (neighbors_of_a.find(node) == neighbors_of_a.end()) {
                        filtered_available_b.push_back(node);
                    }
                }
                for (uint32_t node : all_nodes_b) {
                    if (neighbors_of_a.find(node) == neighbors_of_a.end()) {
                        filtered_all_b.push_back(node);
                    }
                }
                
                // Select second endpoint
                uint32_t selected_b;
                bool b_has_budget = false;
                
                if (!filtered_available_b.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, filtered_available_b.size() - 1);
                    selected_b = filtered_available_b[dist(gen)];
                    b_has_budget = true;
                } else if (!filtered_all_b.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, filtered_all_b.size() - 1);
                    selected_b = filtered_all_b[dist(gen)];
                } else {
                    std::cerr << "[Cluster " << task.cluster_id 
                              << "] Warning: Could not find non-adjacent nodes" << std::endl;
                    break;
                }
                
                uint32_t u = std::min(selected_a, selected_b);
                uint32_t v = std::max(selected_a, selected_b);
                
                uint64_t edge_key = statics::encode_edge(u, v);
                if (planned_edges.find(edge_key) != planned_edges.end()) {
                    continue;
                }
                
                edges_to_add.emplace_back(u, v);
                planned_edges.insert(edge_key);
                current_degrees[u]++;
                current_degrees[v]++;
                connectivity_edges++;
                
                // Consume budgets
                if (has_budget[u]) {
                    uint64_t u_id = g.id_map[u];
                    task.consume_local_degree(u_id, 1);
                    if (task.get_local_available_degree(u_id) == 0) {
                        has_budget[u] = false;
                    }
                    budget_consumed++;
                }
                
                if (has_budget[v]) {
                    uint64_t v_id = g.id_map[v];
                    task.consume_local_degree(v_id, 1);
                    if (task.get_local_available_degree(v_id) == 0) {
                        has_budget[v] = false;
                    }
                    budget_consumed++;
                }
            }
        }
        
        std::cout << "[Cluster " << task.cluster_id << "] Added " << connectivity_edges 
                  << " connectivity edges" << std::endl;
    }
    
    // ==================================================================
    // PHASE 2: MIN DEGREE ENFORCEMENT (Randomized like reference)
    // ==================================================================
    
    // Find nodes needing edges
    std::vector<uint32_t> nodes_needing_edges;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (current_degrees[i] < min_degree) {
            nodes_needing_edges.push_back(i);
        }
    }
    
    if (!nodes_needing_edges.empty()) {
        std::cout << "[Cluster " << task.cluster_id << "] " << nodes_needing_edges.size() 
                  << " nodes need more edges for min degree" << std::endl;
        
        // Process nodes in batches
        const size_t BATCH_SIZE = 1000;
        
        for (size_t batch_start = 0; batch_start < nodes_needing_edges.size(); batch_start += BATCH_SIZE) {
            size_t batch_end = std::min(batch_start + BATCH_SIZE, nodes_needing_edges.size());
            
            for (size_t idx = batch_start; idx < batch_end; ++idx) {
                uint32_t u = nodes_needing_edges[idx];
                
                if (current_degrees[u] >= min_degree) continue;
                
                // Build neighbor exclusion set
                std::vector<bool> is_neighbor(g.num_nodes, false);
                is_neighbor[u] = true;
                
                for (uint32_t j = g.row_ptr[u]; j < g.row_ptr[u + 1]; ++j) {
                    is_neighbor[g.col_idx[j]] = true;
                }
                
                for (const auto& edge : edges_to_add) {
                    if (edge.first == u) is_neighbor[edge.second] = true;
                    if (edge.second == u) is_neighbor[edge.first] = true;
                }
                
                // Add edges until satisfied
                while (current_degrees[u] < min_degree) {
                    if (std::chrono::steady_clock::now() - start_time > MAX_TIME) {
                        std::cout << "[Cluster " << task.cluster_id << "] Timeout during degree enforcement" << std::endl;
                        goto timeout_exit;
                    }
                    
                    uint32_t selected_v = UINT32_MAX;
                    
                    // STRATEGY 1: Try available nodes with budget first
                    std::vector<uint32_t> valid_available;
                    for (uint32_t v : available_nodes) {
                        if (!is_neighbor[v] && has_budget[v]) {
                            valid_available.push_back(v);
                        }
                    }
                    
                    if (!valid_available.empty()) {
                        std::uniform_int_distribution<size_t> dist(0, valid_available.size() - 1);
                        selected_v = valid_available[dist(gen)];
                    } else {
                        // STRATEGY 2: Try any non-neighbor
                        std::vector<uint32_t> all_candidates;
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
                        std::cerr << "[Cluster " << task.cluster_id << "] Warning: Node " << u 
                                  << " cannot reach degree " << min_degree << std::endl;
                        break;
                    }
                    
                    // Add edge
                    uint32_t edge_u = std::min(u, selected_v);
                    uint32_t edge_v = std::max(u, selected_v);
                    
                    edges_to_add.emplace_back(edge_u, edge_v);
                    planned_edges.insert(statics::encode_edge(edge_u, edge_v));
                    
                    current_degrees[u]++;
                    current_degrees[selected_v]++;
                    degree_edges++;
                    
                    // Update budgets
                    if (has_budget[u]) {
                        uint64_t u_id = g.id_map[u];
                        task.consume_local_degree(u_id, 1);
                        if (task.get_local_available_degree(u_id) == 0) {
                            has_budget[u] = false;
                        }
                        budget_consumed++;
                    }
                    
                    if (has_budget[selected_v]) {
                        uint64_t v_id = g.id_map[selected_v];
                        task.consume_local_degree(v_id, 1);
                        if (task.get_local_available_degree(v_id) == 0) {
                            has_budget[selected_v] = false;
                        }
                        budget_consumed++;
                    }
                    
                    is_neighbor[selected_v] = true;
                }
            }
        }
    }
    
timeout_exit:
    
    // Actually add the edges to the graph
    if (!edges_to_add.empty()) {
        add_edges_batch(g, edges_to_add);
    }
    
    std::cout << "[Cluster " << task.cluster_id << "] Budget-aware enforcement complete:" << std::endl;
    std::cout << "  - Total edges added: " << edges_to_add.size() << std::endl;
    std::cout << "  - Connectivity edges: " << connectivity_edges << std::endl;
    std::cout << "  - Degree edges: " << degree_edges << std::endl;
    std::cout << "  - Budget consumed: " << budget_consumed << std::endl;
    
    // Final connectivity check
    auto final_components = find_connected_components_budget(g);
    std::cout << "[Cluster " << task.cluster_id << "] Final: " << final_components.size() 
              << " connected components" << std::endl;
    
    // Count nodes that didn't reach min degree
    uint32_t under_degree_nodes = 0;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (g.get_degree(i) < min_degree) {
            under_degree_nodes++;
        }
    }
    if (under_degree_nodes > 0) {
        std::cout << "[Cluster " << task.cluster_id << "] Warning: " << under_degree_nodes 
                  << " nodes below min degree" << std::endl;
    }
}

#endif // ENFORCE_DEGREE_CONN_WITH_BUDGET_H
