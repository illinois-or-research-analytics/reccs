#ifndef DEG_SEQ_MATCHING_H
#define DEG_SEQ_MATCHING_H

#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

void match_degree_sequence(
    Graph& g, 
    const std::shared_ptr<const std::vector<uint32_t>>& degree_sequence) {
    // Check if the degree sequence is valid
    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes in graph." << std::endl;
        return;
    }

    // Create a degree deficit vector
    // Get current degrees using NodeDegree struct
    std::vector<NodeDegree> node_degrees;
    node_degrees.reserve(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        node_degrees.emplace_back(NodeDegree{i, g.get_degree(i)});
    }

    // Sort by degree
    std::sort(node_degrees.begin(), node_degrees.end(), 
        [](const NodeDegree& a, const NodeDegree& b) {
            return a.degree < b.degree;
        });

    // Create sorted target degree sequence
    std::vector<uint32_t> sorted_target_degrees = *degree_sequence;
    std::sort(sorted_target_degrees.begin(), sorted_target_degrees.end());

    // Compute degree deficit vector using NodeDegree objects and create max heap
    std::vector<NodeDegree> degree_deficit;
    degree_deficit.reserve(g.num_nodes);
    for (uint32_t i = 0; i < g.num_nodes; ++i) {
        int32_t deficit = sorted_target_degrees[i] - node_degrees[i].degree;
        degree_deficit.emplace_back(NodeDegree{node_degrees[i].node, deficit});
    }
    
    // Create max heap based on degree deficit (higher deficit = higher priority)
    std::priority_queue<NodeDegree, std::vector<NodeDegree>, std::less<NodeDegree>> deficit_heap(degree_deficit.begin(), degree_deficit.end());

    // Initialize edge lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);

    // Track edges to add (as pairs where first < second)
    std::set<std::pair<uint32_t, uint32_t>> edges_to_add;

    // Iteratively resolve deficits
    while (!deficit_heap.empty()) {
        // Get node with highest deficit
        NodeDegree current = deficit_heap.top();
        deficit_heap.pop();

        if (current.degree <= 0) {
            // No more deficits to resolve
            break;
        }

        // Find candidates to connect with
        for (uint32_t i = 0; i < g.num_nodes && current.degree > 0; ++i) {
            if (i == current.node || node_degrees[i].degree >= sorted_target_degrees[i]) {
                continue;  // Skip self and already satisfied nodes
            }

            if (!edge_exists_fast(current.node, i)) {
                // Add edge to set
                uint32_t u = std::min(current.node, i);
                uint32_t v = std::max(current.node, i);
                edges_to_add.insert({u, v});
                
                node_degrees[current.node].degree++;
                node_degrees[i].degree++;
                current.degree--;
                
                // Update deficit for the connected node
                deficit_heap.push(NodeDegree{i, sorted_target_degrees[i] - node_degrees[i].degree});
            }
        }
    }

    // Add all edges in batch
    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add_vec(edges_to_add.begin(), edges_to_add.end());
    add_edges_batch(g, edges_to_add_vec);
    std::cout << "Degree sequence matching completed. Added " 
              << edges_to_add.size() << " edges." << std::endl;
}

#endif
