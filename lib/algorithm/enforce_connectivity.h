#ifndef ENFORCE_CONNECTIVITY_H
#define ENFORCE_CONNECTIVITY_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <set>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

// Find connected components using BFS
std::vector<std::vector<uint32_t>> find_connected_components(const Graph& g) {
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

// Get degree of a node
uint32_t get_degree_connectivity(const Graph& g, uint32_t node) {
    return g.row_ptr[node + 1] - g.row_ptr[node];
}

// Enforce connectivity - ensure graph is a single connected component
void enforce_connectivity(Graph& g, uint32_t min_degree) {
    std::cout << "Starting connectivity enforcement on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

    // Find all connected components
    auto components = find_connected_components(g);
    
    // If already connected, we're done
    if (components.size() <= 1) {
        std::cout << "Graph is already connected" << std::endl;
        return;
    }
    
    std::cout << "Found " << components.size() << " components, connecting them..." << std::endl;
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add (as pairs where first < second)
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Track current degrees (will update as we add edges)
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = get_degree_connectivity(g, i);
    }
    
    // Strategy: Build spine through minimum degree nodes from each component
    
    // Create min heap of nodes from each component
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, std::greater<NodeDegree>> min_heap;
    std::vector<uint32_t> component_of_node(g.num_nodes);
    
    for (size_t comp_id = 0; comp_id < components.size(); ++comp_id) {
        // Find minimum degree node in this component
        uint32_t min_node = components[comp_id][0];
        uint32_t min_deg = current_degrees[min_node];
        
        for (uint32_t node : components[comp_id]) {
            component_of_node[node] = comp_id;
            if (current_degrees[node] < min_deg) {
                min_node = node;
                min_deg = current_degrees[node];
            }
        }
        
        min_heap.push({min_node, min_deg});
    }
    
    // Build spine through minimum degree nodes
    std::vector<uint32_t> spine_nodes;
    std::unordered_set<uint32_t> used_components;
    
    while (!min_heap.empty()) {
        NodeDegree nd = min_heap.top();
        min_heap.pop();
        
        uint32_t comp_id = component_of_node[nd.node];
        
        // Skip if we already have a node from this component in the spine
        if (used_components.count(comp_id) > 0) continue;
        
        spine_nodes.push_back(nd.node);
        used_components.insert(comp_id);
    }
    
    // Connect spine nodes to form a single connected component
    for (size_t i = 0; i + 1 < spine_nodes.size(); ++i) {
        uint32_t u = spine_nodes[i];
        uint32_t v = spine_nodes[i + 1];
        
        if (!edge_exists_fast(u, v)) {
            if (u > v) std::swap(u, v);
            edges_to_add.insert({u, v});
            current_degrees[u]++;
            current_degrees[v]++;
        }
    }
    
    // Actually add the edges to the graph using batch addition
    std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_vector);
    
    std::cout << "Added " << edges_to_add.size() << " edges for connectivity" << std::endl;
}

#endif // ENFORCE_CONNECTIVITY_H
