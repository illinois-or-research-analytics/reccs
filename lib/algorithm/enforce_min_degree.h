#ifndef ENFORCE_MIN_DEGREE_H
#define ENFORCE_MIN_DEGREE_H

#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <set>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

// Get degree of a node
uint32_t get_degree(const Graph& g, uint32_t node) {
    return g.row_ptr[node + 1] - g.row_ptr[node];
}

// Enforce minimum degree requirement on all nodes
void enforce_min_degree(Graph& g, uint32_t min_degree) {
    std::cout << "Starting minimum degree enforcement on cluster " 
              << g.id << ". Minimum degree: " << min_degree << std::endl;

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
    
    // Check if all degrees are already satisfied
    bool all_satisfied = true;
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        if (get_degree(g, i) < min_degree) {
            all_satisfied = false;
            break;
        }
    }
    if (all_satisfied) return;
    
    // Build hash set of existing edges for O(1) lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);
    
    // Track edges to add (as pairs where first < second)
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
    
    // Track current degrees (will update as we add edges)
    std::vector<uint32_t> current_degrees(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        current_degrees[i] = get_degree(g, i);
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
    
    // Connect low-degree nodes, prioritizing connections between low-degree nodes
    size_t iterations = 0;
    const size_t MAX_ITERATIONS = g.num_nodes * 10; // Prevent infinite loops
    
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
        
        // If degree already satisfied (due to connections from other nodes), skip
        if (current_degrees[u] >= min_degree) continue;
        
        // Try to connect to other low-degree nodes first
        std::vector<NodeDegree> temp_nodes;
        
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
            }
        }
        
        // Re-add u if still needs edges (shouldn't happen, but just in case)
        if (current_degrees[u] < min_degree) {
            degree_heap.push({u, current_degrees[u]});
            in_heap.insert(u);
        }
    }
    
    // Actually add the edges to the graph using batch addition
    std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_vector);
    
    std::cout << "Added " << edges_to_add.size() << " edges for minimum degree " 
              << min_degree << std::endl;
}

#endif // ENFORCE_MIN_DEGREE_H
