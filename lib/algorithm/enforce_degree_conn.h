#ifndef COMBINED_DEGREE_CONNECTIVITY_H
#define COMBINED_DEGREE_CONNECTIVITY_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <random>
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
uint32_t get_degree(const Graph& g, uint32_t node) {
    return g.row_ptr[node + 1] - g.row_ptr[node];
}

// Combined degree enforcement and connectivity using min heap strategy
void enforce_degree_and_connectivity(Graph& g, uint32_t min_degree) {
    // Check if the requirement is even possible
    if (min_degree >= g.num_nodes) {
        std::cerr << "Error: min_degree " << min_degree 
                  << " >= num_nodes " << g.num_nodes << std::endl;
        return;
    }
    
    // For small graphs with high min_degree, check if it's possible
    if (g.num_nodes > 1 && min_degree > g.num_nodes - 1) {
        std::cerr << "Error: Impossible to achieve min_degree " << min_degree 
                  << " with " << g.num_nodes << " nodes" << std::endl;
        return;
    }
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add (as pairs where first < second)
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // First, find all connected components
    auto components = find_connected_components(g);
    
    // If already connected and all degrees satisfied, we're done
    if (components.size() == 1) {
        bool all_satisfied = true;
        for (uint32_t i = 0; i < g.num_nodes; ++i) {
            if (get_degree(g, i) < min_degree) {
                all_satisfied = false;
                break;
            }
        }
        if (all_satisfied) return;
    }
    
    // Track current degrees (will update as we add edges)
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = get_degree(g, i);
    }
    
    // Strategy: Build spine through minimum degree nodes
    
    // Step 1: Create min heap of nodes from each component
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
    
    // Step 2: Build spine through minimum degree nodes
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
    
    // Connect spine nodes
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
    
    // Step 3: Build heap of all nodes that still need more edges
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, std::greater<NodeDegree>> degree_heap;
    std::unordered_set<uint32_t> in_heap;
    
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (current_degrees[i] < min_degree) {
            degree_heap.push({i, current_degrees[i]});
            in_heap.insert(i);
        }
    }
    
    // Step 4: Connect low-degree nodes, prioritizing connections between low-degree nodes
    size_t iterations = 0;
    const size_t MAX_ITERATIONS = g.num_nodes * 10; // Prevent infinite loops
    
    while (!degree_heap.empty() && iterations < MAX_ITERATIONS) {
        iterations++;
        if (iterations % 1000 == 0) {
            std::cout << "Progress: " << iterations << " iterations, " 
                      << degree_heap.size() << " nodes remaining" << std::endl;
        }
        NodeDegree nd = degree_heap.top();
        degree_heap.pop();
        uint32_t u = nd.node;
        in_heap.erase(u);
        
        // If degree already satisfied (due to connections from other nodes), skip
        if (current_degrees[u] >= min_degree) continue;
        
        // Try to connect to other low-degree nodes first
        std::vector<NodeDegree> temp_nodes;
        bool found_connections = false;
        
        while (!degree_heap.empty() && current_degrees[u] < min_degree) {
            NodeDegree other = degree_heap.top();
            degree_heap.pop();
            uint32_t v = other.node;
            
            if (!edge_exists_fast(u, v) && 
                edges_to_add.count({std::min(u,v), std::max(u,v)}) == 0) {
                
                // Add edge
                edges_to_add.insert({std::min(u,v), std::max(u,v)});
                current_degrees[u]++;
                current_degrees[v]++;
                found_connections = true;
                
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
        
        // If still need edges, connect to any available nodes
        if (current_degrees[u] < min_degree) {
            std::vector<uint32_t> candidates;
            
            for (uint32_t v = 0; v < g.num_nodes; ++v) {
                if (v != u && !edge_exists_fast(u, v) && 
                    edges_to_add.count({std::min(u,v), std::max(u,v)}) == 0) {
                    candidates.push_back(v);
                }
            }
            
            // Check if we have enough candidates
            if (candidates.size() < min_degree - current_degrees[u]) {
                std::cerr << "Warning: Node " << u << " cannot reach degree " << min_degree
                          << " (only " << candidates.size() << " candidates available)" << std::endl;
            }
            
            // Sort candidates by degree (prefer lower degree nodes)
            std::sort(candidates.begin(), candidates.end(), 
                     [&current_degrees](uint32_t a, uint32_t b) {
                         return current_degrees[a] < current_degrees[b];
                     });
            
            // Add edges until degree satisfied
            for (uint32_t v : candidates) {
                if (current_degrees[u] >= min_degree) break;
                
                edges_to_add.insert({std::min(u,v), std::max(u,v)});
                current_degrees[u]++;
                current_degrees[v]++;
                
                // If v was in heap and now satisfied, we'll skip it when popped
            }
        }
        
        // Re-add u if still needs edges (shouldn't happen, but just in case)
        if (current_degrees[u] < min_degree) {
            degree_heap.push({u, current_degrees[u]});
            in_heap.insert(u);
        }
    }
    
    // Step 5: Actually add the edges to the graph using batch addition
    std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_vector);
    
    std::cout << "Added " << edges_to_add.size() << " edges for degree " 
              << min_degree << " and connectivity" << std::endl;

    // DEBUG CODE
    // save_graph_edgelist(g.id + "_min_degree_" + std::to_string(min_degree) + ".tsv", g, true);
    // END DEBUG CODE
}

#endif // COMBINED_DEGREE_CONNECTIVITY_H
