#ifndef DEG_SEQ_MATCHING_H
#define DEG_SEQ_MATCHING_H

#include <vector>
#include <memory>
#include <iostream>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"

void match_degree_sequence(
    Graph& g, 
    const std::shared_ptr<const std::vector<uint32_t>>& degree_sequence) {
    // Verify that degree_sequence is sorted in non-decreasing order
    if (!degree_sequence || degree_sequence->empty()) {
        std::cerr << "Error: Degree sequence is null or empty." << std::endl;
        return;
    }

    // Check if the degree sequence is valid
    if (degree_sequence->size() != g.num_nodes) {
        std::cerr << "Error: Degree sequence size does not match number of nodes in graph." 
                  << "Expected: " << g.num_nodes
                  << ", Got: " << degree_sequence->size() << std::endl;
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
    std::make_heap(degree_deficit.begin(), degree_deficit.end(), 
        [](const NodeDegree& a, const NodeDegree& b) {
            return a.degree < b.degree;  // For max heap, use less-than comparator
        });

    // Initialize edge lookup
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists_fast = statics::create_edge_exists_checker(existing_edges);

    // TODO: Implement the matching algorithm
    while (!degree_deficit.empty()) {
        // Get node with highest deficit
        std::pop_heap(degree_deficit.begin(), degree_deficit.end(), 
            [](const NodeDegree& a, const NodeDegree& b) {
                return a.degree < b.degree;  // For max heap, use less-than comparator
            });
        NodeDegree current = degree_deficit.back();
        degree_deficit.pop_back();

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
                // Add edge
                g.add_edge(current.node, i);
                node_degrees[current.node].degree++;
                node_degrees[i].degree++;
                current.degree--;
                
                // Update deficit for the connected node
                degree_deficit.push_back(NodeDegree{i, sorted_target_degrees[i] - node_degrees[i].degree});
                std::push_heap(degree_deficit.begin(), degree_deficit.end(), 
                    [](const NodeDegree& a, const NodeDegree& b) {
                        return a.degree < b.degree;  // For max heap, use less-than comparator
                    });
            }
        }
    }
}

#endif
