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

uint32_t get_degree_budget(const Graph& g, uint32_t node) {
    return g.row_ptr[node + 1] - g.row_ptr[node];
}

/**
 * Degree-aware minimum degree enforcement
 * Replicates Python logic: prioritize available degree nodes when adding edges
 */
void enforce_min_degree_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;
    
    std::cout << "Starting minimum degree enforcement with budget on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    if (min_degree >= g.num_nodes) {
        std::cerr << "Error: min_degree " << min_degree 
                  << " >= num_nodes " << g.num_nodes << std::endl;
        return;
    }
    
    if (g.num_nodes > 1 && min_degree > g.num_nodes - 1) {
        std::cerr << "Error: Impossible to achieve min_degree " << min_degree 
                  << " with " << g.num_nodes << " nodes" << std::endl;
        return;
    }
    
    // Check if all degrees are already satisfied
    bool all_satisfied = true;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (get_degree_budget(g, i) < min_degree) {
            all_satisfied = false;
            break;
        }
    }
    if (all_satisfied) return;
    
    // Get available nodes for this cluster (intersection with cluster nodes)
    auto cluster_available_nodes = degree_aware_stages::get_cluster_available_nodes(task);
    
    std::cout << "Cluster has " << cluster_available_nodes.size() 
              << " nodes with available degree budget" << std::endl;
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Track current degrees (will update as we add edges)
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = get_degree_budget(g, i);
    }
    
    // Build heap of all nodes that need more edges
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, std::greater<NodeDegree>> degree_heap;
    std::unordered_set<uint32_t> in_heap;
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (current_degrees[i] < min_degree) {
            degree_heap.push({i, current_degrees[i]});
            in_heap.insert(i);
        }
    }
    
    // Random number generator for tie-breaking
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    size_t iterations = 0;
    const size_t MAX_ITERATIONS = g.num_nodes * 10;
    
    while (!degree_heap.empty() && iterations < MAX_ITERATIONS) {
        iterations++;
        if (iterations % 1000 == 0) {
            std::cout << "Min degree progress: " << iterations << " iterations, " 
                      << degree_heap.size() << " nodes remaining" << std::endl;
        }
        
        NodeDegree nd = degree_heap.top();
        degree_heap.pop();
        uint32_t u = nd.node;
        in_heap.erase(u);
        
        // If degree already satisfied, skip
        if (current_degrees[u] >= min_degree) continue;
        
        // Get current neighbors of u
        std::unordered_set<uint32_t> neighbors;
        for (uint32_t j = g.row_ptr[u]; j < g.row_ptr[u + 1]; ++j) {
            neighbors.insert(g.col_idx[j]);
        }
        neighbors.insert(u); // Self-loop prevention
        
        // Add edges from planned additions
        for (const auto& edge : edges_to_add) {
            if (edge.first == u) neighbors.insert(edge.second);
            if (edge.second == u) neighbors.insert(edge.first);
        }
        
        // Build candidate lists: available vs regular nodes
        std::vector<uint32_t> available_candidates;
        std::vector<uint32_t> regular_candidates;
        
        for (uint32_t v = 0; v < g.num_nodes; ++v) {
            if (neighbors.count(v) > 0) continue; // Already connected or self
            
            uint64_t v_id = g.id_map[v];
            if (cluster_available_nodes.count(v_id) > 0 && 
                task.degree_manager->get_available_degree(v_id) > 0) {
                available_candidates.push_back(v);
            } else {
                regular_candidates.push_back(v);
            }
        }
        
        // Try to connect to other low-degree nodes first, prioritizing available degree nodes
        std::vector<NodeDegree> temp_nodes;
        
        while (!degree_heap.empty() && current_degrees[u] < min_degree) {
            NodeDegree other = degree_heap.top();
            degree_heap.pop();
            uint32_t v = other.node;
            
            if (neighbors.count(v) == 0) { // Can connect
                bool v_has_budget = false;
                uint64_t v_id = g.id_map[v];
                if (cluster_available_nodes.count(v_id) > 0) {
                    v_has_budget = task.degree_manager->get_available_degree(v_id) > 0;
                }
                
                // Add edge
                edges_to_add.insert({std::min(u,v), std::max(u,v)});
                current_degrees[u]++;
                current_degrees[v]++;
                total_edges_added++;
                
                // Track if this edge uses available degree budget
                uint64_t u_id = g.id_map[u];
                bool u_has_budget = cluster_available_nodes.count(u_id) > 0 && 
                                   task.degree_manager->get_available_degree(u_id) > 0;
                
                if (u_has_budget || v_has_budget) {
                    degree_corrected_edges++;
                }
                
                // Update neighbors set
                neighbors.insert(v);
                
                // Re-add v to heap if it still needs edges
                if (current_degrees[v] < min_degree) {
                    temp_nodes.push_back({v, current_degrees[v]});
                } else {
                    in_heap.erase(v);
                }
            } else {
                // Can't connect to this node, save it
                temp_nodes.push_back(other);
            }
        }
        
        // Re-add saved nodes to heap
        for (const auto& node : temp_nodes) {
            degree_heap.push(node);
        }
        
        // If still need edges, connect to available candidates first, then regular
        while (current_degrees[u] < min_degree) {
            uint32_t v = UINT32_MAX;
            bool found_candidate = false;
            
            // Try available candidates first
            if (!available_candidates.empty()) {
                v = degree_aware_stages::select_edge_endpoint(available_candidates, g, task, gen);
                if (v != UINT32_MAX) {
                    // Remove from available candidates
                    available_candidates.erase(
                        std::remove(available_candidates.begin(), available_candidates.end(), v),
                        available_candidates.end());
                    found_candidate = true;
                }
            }
            
            // If no available candidates, try regular candidates
            if (!found_candidate && !regular_candidates.empty()) {
                v = degree_aware_stages::select_edge_endpoint(regular_candidates, g, task, gen);
                if (v != UINT32_MAX) {
                    // Remove from regular candidates
                    regular_candidates.erase(
                        std::remove(regular_candidates.begin(), regular_candidates.end(), v),
                        regular_candidates.end());
                    found_candidate = true;
                }
            }
            
            if (!found_candidate) {
                std::cerr << "Warning: Node " << u << " cannot reach degree " << min_degree
                          << " (no more candidates available)" << std::endl;
                break;
            }
            
            // Add the edge
            edges_to_add.insert({std::min(u,v), std::max(u,v)});
            current_degrees[u]++;
            current_degrees[v]++;
            total_edges_added++;
            
            // Track if this edge uses available degree budget
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[v];
            bool u_has_budget = cluster_available_nodes.count(u_id) > 0 && 
                               task.degree_manager->get_available_degree(u_id) > 0;
            bool v_has_budget = cluster_available_nodes.count(v_id) > 0 && 
                               task.degree_manager->get_available_degree(v_id) > 0;
            
            if (u_has_budget || v_has_budget) {
                degree_corrected_edges++;
            }
            
            neighbors.insert(v);
        }
        
        // Re-add u if still needs edges (shouldn't happen, but safety check)
        if (current_degrees[u] < min_degree) {
            degree_heap.push({u, current_degrees[u]});
            in_heap.insert(u);
        }
    }
    
    // Actually add the edges to the graph using batch addition
    std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_vector);
    
    // Consume degree budgets for edges that were added
    for (const auto& edge : edges_vector) {
        uint64_t u_id = g.id_map[edge.first];
        uint64_t v_id = g.id_map[edge.second];
        
        // Try to consume budget (thread-safe)
        task.degree_manager->try_consume_degree(u_id, 1);
        task.degree_manager->try_consume_degree(v_id, 1);
    }
    
    std::cout << "Added " << edges_to_add.size() << " edges for minimum degree " 
              << min_degree << ". Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%)" << std::endl;
}

#endif // ENFORCE_MIN_DEGREE_WITH_BUDGET_H
