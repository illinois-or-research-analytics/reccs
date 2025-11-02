#ifndef DEG_SEQ_MATCHING_PP_H
#define DEG_SEQ_MATCHING_PP_H

#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>

#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"
#include "../data_structures/available_node_degrees.h"

/**
 * Matches the degree sequence of a graph using pre-computed degree deficits.
 * Heavily optimized version that eliminates unnecessary copying.
 */
void match_degree_sequence_pp(GraphTask& task) {
    // Extract graph and degree manager from task
    Graph& g = *task.subgraph;
    auto degree_manager = task.degree_manager;

    if (!degree_manager) {
        std::cerr << "Error: Task's degree manager is null." << std::endl;
        return;
    }

    // Get available nodes for this cluster
    auto available_nodes = task.get_local_available_nodes();

    if (available_nodes.empty()) {
        std::cout << "[Cluster " << task.cluster_id << "]: No nodes need additional edges." << std::endl;
        return;
    }

    // Python: available_node_degrees dict and available_node_set
    std::unordered_map<uint32_t, uint32_t> available_node_degrees; // local_idx -> deficit

    // OPTIMIZATION: Use bitmap instead of unordered_set - can be memcpy'd!
    std::vector<uint8_t> is_available_bitmap(g.num_nodes, 0);
    std::vector<uint32_t> available_node_list; // For iteration

    // Initialize from degree deficits
    for (uint64_t global_node_id : available_nodes) {
        auto it = g.node_map.find(global_node_id);
        if (it != g.node_map.end()) {
            uint32_t local_idx = it->second;
            int32_t deficit = task.get_local_available_degree(global_node_id);
            if (deficit > 0) {
                available_node_degrees[local_idx] = static_cast<uint32_t>(deficit);
                is_available_bitmap[local_idx] = 1;
                available_node_list.push_back(local_idx);
            }
        }
    }

    if (available_node_degrees.empty()) {
        std::cout << "[Cluster " << task.cluster_id << "]: No nodes in this cluster need additional edges." << std::endl;
        return;
    }

    // Python: exist_neighbor dict
    // Build persistent neighbor map ONCE at start (NOT included in Python's main loop timing!)
    auto setup_start = std::chrono::high_resolution_clock::now();
    std::cout << "[Cluster " << task.cluster_id << "]: Building neighbor map..." << std::endl;

    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> exist_neighbor;
    exist_neighbor.reserve(g.num_nodes);

    for (uint32_t node = 0; node < g.num_nodes; ++node) {
        exist_neighbor[node] = std::unordered_set<uint32_t>();
        for (uint32_t idx = g.row_ptr[node]; idx < g.row_ptr[node + 1]; ++idx) {
            exist_neighbor[node].insert(g.col_idx[idx]);
        }
    }

    auto setup_end = std::chrono::high_resolution_clock::now();
    auto setup_duration = std::chrono::duration_cast<std::chrono::milliseconds>(setup_end - setup_start);
    std::cout << "[Cluster " << task.cluster_id << "]: Neighbor map build time: "
              << setup_duration.count() / 1000.0 << " seconds (NOT in main loop)" << std::endl;

    // Python: degree_edges = set() - Collect edges to add in batch
    std::vector<std::pair<uint32_t, uint32_t>> degree_edges;
    degree_edges.reserve(available_node_list.size() * 10); // Rough estimate

    // Python: max_heap using heapq (max heap with negative degrees)
    auto heap_comparator = [](const NodeDegree& a, const NodeDegree& b) {
        if (a.degree != b.degree) return a.degree < b.degree; // Max heap
        return a.node > b.node;
    };

    std::priority_queue<NodeDegree, std::vector<NodeDegree>, decltype(heap_comparator)>
        max_heap(heap_comparator);

    for (const auto& [local_idx, degree] : available_node_degrees) {
        max_heap.emplace(NodeDegree{local_idx, static_cast<int32_t>(degree)});
    }

    uint32_t nodes_processed = 0;
    std::cout << "[Cluster " << task.cluster_id << "]: Processing " << max_heap.size() << " nodes..." << std::endl;

    // CRITICAL OPTIMIZATION: Use bitmaps that can be fast-copied (memcpy-able!)
    std::vector<uint8_t> available_non_neighbors_bitmap(g.num_nodes);
    std::vector<uint32_t> available_non_neighbors_vec;
    available_non_neighbors_vec.reserve(available_node_list.size());

    // START TIMING: Main degree matching loop
    auto loop_start = std::chrono::high_resolution_clock::now();

    // Main loop - optimized version
    while (!max_heap.empty()) {
        // Line 271: heapq.heappop(max_heap)
        NodeDegree current = max_heap.top();
        max_heap.pop();

        uint32_t available_c_node = current.node;

        // Line 274-275: Check if node is still available
        if (available_node_degrees.find(available_c_node) == available_node_degrees.end()) {
            continue;
        }

        // Line 278-284: MATCH PYTHON with fast bitmap copy (like memcpy!)
        // Python: available_non_neighbors = available_node_set.copy()
        //         for neighbor in neighbors: available_non_neighbors.discard(neighbor)

        auto neighbor_it = exist_neighbor.find(available_c_node);
        if (neighbor_it == exist_neighbor.end()) {
            available_node_degrees.erase(available_c_node);
            is_available_bitmap[available_c_node] = 0;
            continue;
        }

        auto& neighbor_set = neighbor_it->second;

        // PYTHON'S EXACT ALGORITHM with fast bitmap:
        // 1. Copy bitmap (optimized memcpy in vector copy)
        available_non_neighbors_bitmap = is_available_bitmap;

        // 2. Erase self and neighbors
        available_non_neighbors_bitmap[available_c_node] = 0;
        for (uint32_t neighbor : neighbor_set) {
            available_non_neighbors_bitmap[neighbor] = 0;
        }

        // 3. Build vector of available non-neighbors for pop operations
        available_non_neighbors_vec.clear();
        for (uint32_t node : available_node_list) {
            if (available_non_neighbors_bitmap[node]) {
                available_non_neighbors_vec.push_back(node);
            }
        }

        // Line 287-290: Compute avail_k
        uint32_t avail_k = std::min(
            available_node_degrees[available_c_node],
            static_cast<uint32_t>(available_non_neighbors_vec.size())
        );

        // Line 293-309: Add edges
        for (uint32_t i = 0; i < avail_k; ++i) {
            if (available_non_neighbors_vec.empty()) break;

            // Line 295: available_non_neighbors.pop() - O(1)
            uint32_t edge_end = available_non_neighbors_vec.back();
            available_non_neighbors_vec.pop_back();

            // Line 298: degree_edges.add((available_c_node, edge_end))
            degree_edges.emplace_back(available_c_node, edge_end);

            // Line 301-302: Update exist_neighbor
            // We already have neighbor_it from above
            neighbor_it->second.insert(edge_end);
            exist_neighbor[edge_end].insert(available_c_node);

            // Line 305-308: Update degree of edge_end
            available_node_degrees[edge_end]--;
            if (available_node_degrees[edge_end] == 0) {
                is_available_bitmap[edge_end] = 0;
                available_node_degrees.erase(edge_end);
            }
        }

        // Line 311-312: Remove processed node
        available_node_degrees.erase(available_c_node);
        is_available_bitmap[available_c_node] = 0;

        // Line 313-318: Progress logging
        nodes_processed++;
        if (nodes_processed % 1000 == 0) {
            std::cout << "[Cluster " << task.cluster_id << "] Processed " << nodes_processed << " nodes" << std::endl;
        }
    }

    // END TIMING: Main degree matching loop
    auto loop_end = std::chrono::high_resolution_clock::now();
    auto loop_duration = std::chrono::duration_cast<std::chrono::milliseconds>(loop_end - loop_start);

    std::cout << "[Cluster " << task.cluster_id << "]: Main loop time: "
              << loop_duration.count() / 1000.0 << " seconds" << std::endl;

    // Add all edges in one efficient batch operation
    auto batch_start = std::chrono::high_resolution_clock::now();
    std::cout << "[Cluster " << task.cluster_id << "]: Adding " << degree_edges.size() << " edges to graph..." << std::endl;
    add_edges_batch(g, degree_edges);
    auto batch_end = std::chrono::high_resolution_clock::now();
    auto batch_duration = std::chrono::duration_cast<std::chrono::milliseconds>(batch_end - batch_start);

    std::cout << "[Cluster " << task.cluster_id << "]: Batch add time: "
              << batch_duration.count() / 1000.0 << " seconds" << std::endl;

    std::cout << "[Cluster " << task.cluster_id << "]: TOTAL DEGREE MATCHING TIME: "
              << (loop_duration.count() + batch_duration.count()) / 1000.0 << " seconds" << std::endl;
    std::cout << "[Cluster " << task.cluster_id << "]: Complete. "
              << "Edges added: " << degree_edges.size()
              << ", Nodes processed: " << nodes_processed << std::endl;
}

#endif // DEG_SEQ_MATCHING_PP_H
