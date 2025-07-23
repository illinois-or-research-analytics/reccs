#ifndef ENFORCE_CONNECTIVITY_WITH_BUDGET_H
#define ENFORCE_CONNECTIVITY_WITH_BUDGET_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <set>
#include <random>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"
#include "../data_structures/available_node_degrees.h"

// Find connected components using BFS
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
 * Degree-aware connectivity enforcement using local degree management
 * Each task operates on its own local degree budgets - NO SHARED STATE!
 */
void enforce_connectivity_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;

    auto start_time = std::chrono::steady_clock::now();
    const auto MAX_TIME = std::chrono::seconds(30); // 30 second timeout
    
    std::cout << "Starting local connectivity enforcement on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    // Find all connected components
    auto components = find_connected_components_budget(g);
    
    // If already connected, we're done
    if (components.size() <= 1) {
        std::cout << "Graph is already connected" << std::endl;
        return;
    }
    
    std::cout << "Found " << components.size() << " components, connecting with " 
              << min_degree << " edges..." << std::endl;    

    // // Initialize local degrees for this task (lazy initialization)
    // task.initialize_local_degrees();
    
    // Get local available nodes (no shared state!)
    const auto& local_available_nodes = task.get_local_available_nodes();
    
    std::cout << "Cluster has " << local_available_nodes.size() 
              << " nodes with local available degree budget" << std::endl;

    // Create lookup set for available nodes from LOCAL data
    std::unordered_set<uint64_t> available_node_set(local_available_nodes.begin(),
                                                   local_available_nodes.end());
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add
    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    // Strategy: Draw min_degree_requirement edges between components using random selection
    
    // Pick largest component as the main partitions to connect
    std::sort(components.begin(), components.end(), 
              [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                  return a.size() > b.size();
              });

    // Connect all components to the largest component first (ensure connectivity)
    for (size_t i = 1; i < components.size(); ++i) {
        const auto& component_a = components[0];  // Largest component
        const auto& component_b = components[i];  // Current component
        
        // Draw min_degree_requirement edges between these two components
        for (uint32_t edge_count = 0; edge_count < min_degree; ++edge_count) {
            if (std::chrono::steady_clock::now() - start_time > MAX_TIME) {
                std::cout << "Timeout reached, stopping connectivity enforcement" << std::endl;
                break;
            }
            
            // Build candidate lists for both sides
            std::vector<uint32_t> available_nodes_a, all_nodes_a;
            std::vector<uint32_t> available_nodes_b, all_nodes_b;
            
            // Component A candidates
            for (uint32_t node : component_a) {
                all_nodes_a.push_back(node);
                uint64_t node_id = g.id_map[node];
                if (task.get_local_available_degree(node_id) > 0) {
                    available_nodes_a.push_back(node);
                }
            }
            
            // Component B candidates  
            for (uint32_t node : component_b) {
                all_nodes_b.push_back(node);
                uint64_t node_id = g.id_map[node];
                if (task.get_local_available_degree(node_id) > 0) {
                    available_nodes_b.push_back(node);
                }
            }
            
            // Step 1: Select first endpoint (from component A)
            uint32_t selected_a;
            bool a_has_budget = false;
            
            if (!available_nodes_a.empty()) {
                std::uniform_int_distribution<size_t> dist_a(0, available_nodes_a.size() - 1);
                selected_a = available_nodes_a[dist_a(gen)];
                a_has_budget = true;
            } else {
                std::uniform_int_distribution<size_t> dist_a(0, all_nodes_a.size() - 1);
                selected_a = all_nodes_a[dist_a(gen)];
            }
            
            // Step 2: Apply neighbor exclusion - remove neighbors of selected_a from component B candidates
            std::unordered_set<uint32_t> neighbors_of_a;
            for (uint32_t idx = g.row_ptr[selected_a]; idx < g.row_ptr[selected_a + 1]; ++idx) {
                neighbors_of_a.insert(g.col_idx[idx]);
            }
            
            // Filter available nodes in component B (exclude neighbors)
            std::vector<uint32_t> filtered_available_b;
            for (uint32_t node : available_nodes_b) {
                if (neighbors_of_a.find(node) == neighbors_of_a.end()) {
                    filtered_available_b.push_back(node);
                }
            }
            
            // Filter all nodes in component B (exclude neighbors)  
            std::vector<uint32_t> filtered_all_b;
            for (uint32_t node : all_nodes_b) {
                if (neighbors_of_a.find(node) == neighbors_of_a.end()) {
                    filtered_all_b.push_back(node);
                }
            }
            
            // Step 3: Select second endpoint (from filtered component B)
            uint32_t selected_b;
            bool b_has_budget = false;
            
            if (!filtered_available_b.empty()) {
                std::uniform_int_distribution<size_t> dist_b(0, filtered_available_b.size() - 1);
                selected_b = filtered_available_b[dist_b(gen)];
                b_has_budget = true;
            } else if (!filtered_all_b.empty()) {
                std::uniform_int_distribution<size_t> dist_b(0, filtered_all_b.size() - 1);
                selected_b = filtered_all_b[dist_b(gen)];
            } else {
                std::cerr << "Warning: Could not find non-adjacent nodes to connect components" << std::endl;
                break;
            }
            
            // No need to check if edge exists since neighbor exclusion guarantees it doesn't
            uint32_t u = std::min(selected_a, selected_b);
            uint32_t v = std::max(selected_a, selected_b);
            
            // Check if already in our planned edges
            bool already_planned = false;
            for (const auto& edge : edges_to_add) {
                if (edge.first == u && edge.second == v) {
                    already_planned = true;
                    break;
                }
            }
            
            if (already_planned) {
                continue; // Skip already planned edge
            }
            
            // Add the edge
            edges_to_add.emplace_back(u, v);
            total_edges_added++;
            
            // Track if this edge uses available degree budget using LOCAL data
            if (a_has_budget || b_has_budget) {
                degree_corrected_edges++;
            }
            
            // Consume local budgets (no contention!)
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[v];
            task.consume_local_degree(u_id, 1);
            task.consume_local_degree(v_id, 1);
            
            std::cout << "Connecting component " << 0 << " (node " << selected_a 
                      << ") to component " << i << " (node " << selected_b 
                      << ") - edge " << (edge_count + 1) << "/" << min_degree << std::endl;
        }
    }
    
    // Batch add all edges at once
    if (!edges_to_add.empty()) {
        add_edges_batch(g, edges_to_add);
    }
    
    // Report final local available nodes
    const auto& final_available = task.get_local_available_nodes();
    
    std::cout << "Added " << edges_to_add.size() << " edges for connectivity (" 
              << min_degree << " per component pair). "
              << "Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%). Remaining local budget nodes: " << final_available.size() << std::endl;
}

#endif // ENFORCE_CONNECTIVITY_WITH_BUDGET_H
