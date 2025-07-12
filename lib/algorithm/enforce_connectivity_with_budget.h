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

uint32_t get_degree_connectivity_budget(const Graph& g, uint32_t node) {
    return g.row_ptr[node + 1] - g.row_ptr[node];
}

/**
 * Degree-aware connectivity enforcement
 * Replicates Python logic for Stage 2: connecting disconnected components while 
 * prioritizing nodes with available degree budget
 */
void enforce_connectivity_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_degree = task.min_degree_requirement;
    
    std::cout << "Starting connectivity enforcement with budget on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    // Find all connected components
    auto components = find_connected_components_budget(g);
    
    // If already connected, we're done
    if (components.size() <= 1) {
        std::cout << "Graph is already connected" << std::endl;
        return;
    }
    
    std::cout << "Found " << components.size() << " components, connecting them..." << std::endl;
    
    // Get available nodes for this cluster
    auto cluster_available_nodes = degree_aware_stages::get_cluster_available_nodes(task);
    
    std::cout << "Cluster has " << cluster_available_nodes.size() 
              << " nodes with available degree budget" << std::endl;
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Track current degrees
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = get_degree_connectivity_budget(g, i);
    }
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    // Strategy: Connect components through lowest-degree nodes, preferring available budget nodes
    
    // Create priority queue of representative nodes from each component
    // Priority: available budget nodes with lower degrees first
    struct ComponentRepresentative {
        uint32_t node;
        uint32_t component_id;
        uint32_t degree;
        bool has_budget;
        
        // Comparator: prefer available budget nodes, then lower degree
        bool operator>(const ComponentRepresentative& other) const {
            if (has_budget != other.has_budget) {
                return !has_budget; // has_budget nodes have higher priority (lower in max-heap)
            }
            return degree > other.degree; // Lower degree has higher priority
        }
    };
    
    std::priority_queue<ComponentRepresentative, 
                       std::vector<ComponentRepresentative>, 
                       std::greater<ComponentRepresentative>> comp_queue;
    
    std::vector<uint32_t> component_of_node(g.num_nodes);
    
    // Find representative node for each component (prefer available budget, then min degree)
    for (size_t comp_id = 0; comp_id < components.size(); ++comp_id) {
        uint32_t best_node = components[comp_id][0];
        uint32_t best_degree = current_degrees[best_node];
        bool best_has_budget = false;
        
        // Check if best node has budget
        uint64_t best_node_id = g.id_map[best_node];
        if (cluster_available_nodes.count(best_node_id) > 0) {
            best_has_budget = task.degree_manager->get_available_degree(best_node_id) > 0;
        }
        
        for (uint32_t node : components[comp_id]) {
            component_of_node[node] = comp_id;
            
            uint64_t node_id = g.id_map[node];
            bool has_budget = cluster_available_nodes.count(node_id) > 0 && 
                             task.degree_manager->get_available_degree(node_id) > 0;
            uint32_t node_degree = current_degrees[node];
            
            // Prefer nodes with budget, then lower degree
            bool is_better = false;
            if (has_budget && !best_has_budget) {
                is_better = true;
            } else if (has_budget == best_has_budget) {
                is_better = (node_degree < best_degree);
            }
            
            if (is_better) {
                best_node = node;
                best_degree = node_degree;
                best_has_budget = has_budget;
            }
        }
        
        comp_queue.push({best_node, static_cast<uint32_t>(comp_id), best_degree, best_has_budget});
    }
    
    // Build spine connecting components
    std::vector<uint32_t> spine_nodes;
    std::unordered_set<uint32_t> used_components;
    
    while (!comp_queue.empty()) {
        ComponentRepresentative repr = comp_queue.top();
        comp_queue.pop();
        
        // Skip if we already have a node from this component in the spine
        if (used_components.count(repr.component_id) > 0) continue;
        
        spine_nodes.push_back(repr.node);
        used_components.insert(repr.component_id);
    }
    
    std::cout << "Selected spine nodes: ";
    for (uint32_t node : spine_nodes) {
        uint64_t node_id = g.id_map[node];
        bool has_budget = cluster_available_nodes.count(node_id) > 0 && 
                         task.degree_manager->get_available_degree(node_id) > 0;
        std::cout << node << (has_budget ? "*" : "") << " ";
    }
    std::cout << std::endl;
    
    // Connect spine nodes to form a single connected component
    for (size_t i = 0; i + 1 < spine_nodes.size(); ++i) {
        uint32_t u = spine_nodes[i];
        uint32_t v = spine_nodes[i + 1];
        
        if (!edge_exists_fast(u, v)) {
            if (u > v) std::swap(u, v);
            edges_to_add.insert({u, v});
            current_degrees[u]++;
            current_degrees[v]++;
            total_edges_added++;
            
            // Check if this edge uses available degree budget
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[v];
            bool u_has_budget = cluster_available_nodes.count(u_id) > 0 && 
                               task.degree_manager->get_available_degree(u_id) > 0;
            bool v_has_budget = cluster_available_nodes.count(v_id) > 0 && 
                               task.degree_manager->get_available_degree(v_id) > 0;
            
            if (u_has_budget || v_has_budget) {
                degree_corrected_edges++;
            }
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
    
    std::cout << "Added " << edges_to_add.size() << " edges for connectivity. "
              << "Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%)" << std::endl;
}

#endif // ENFORCE_CONNECTIVITY_WITH_BUDGET_H
