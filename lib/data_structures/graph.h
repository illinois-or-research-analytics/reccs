#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <omp.h>

// CSR representation for undirected graph
struct Graph {
    std::vector<uint32_t> row_ptr;  // Offsets for each node's edge list
    std::vector<uint32_t> col_idx;  // Target nodes
    
    // Node ID mapping (if needed)
    std::unordered_map<uint64_t, uint32_t> node_map;
    std::vector<uint64_t> id_map;
    
    // Graph info
    size_t num_nodes = 0;
    size_t num_edges = 0; // This counts each undirected edge once

    // Add an edge to the graph
    void add_edge(uint32_t from, uint32_t to) {
        // Check if nodes exist
        if (from >= num_nodes || to >= num_nodes) {
            return;
        }
        
        // Check if edge already exists
        for (uint32_t i = row_ptr[from]; i < row_ptr[from + 1]; ++i) {
            if (col_idx[i] == to) {
                return; // Edge already exists
            }
        }
        
        // We need to shift all row pointers after 'from' to accommodate the new edge
        for (uint32_t i = from + 1; i <= num_nodes; ++i) {
            row_ptr[i]++;
        }
        
        // Insert the new edge
        col_idx.insert(col_idx.begin() + row_ptr[from], to);
        
        // Add the reverse edge (for undirected graph)
        for (uint32_t i = to + 1; i <= num_nodes; ++i) {
            row_ptr[i]++;
        }
        
        col_idx.insert(col_idx.begin() + row_ptr[to], from);
        
        // Update edge count
        num_edges++;
    }
};

// Sort adjacency lists in parallel
void sort_adjacency_lists_parallel(Graph& graph, int num_threads, bool verbose = false) {
    if (verbose) {
        std::cout << "Sorting adjacency lists..." << std::endl;
    }
    
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        
        if (end > start) {
            std::sort(&graph.col_idx[start], &graph.col_idx[end]);
        }
    }
}

// Remove self-loops and duplicate edges in parallel
void clean_graph_parallel(Graph& graph, int num_threads, bool verbose = false) {
    if (verbose) {
        std::cout << "Removing self-loops and duplicate edges..." << std::endl;
    }
    
    // First, sort adjacency lists to place duplicates adjacent to each other
    sort_adjacency_lists_parallel(graph, num_threads, verbose);
    
    // Count unique edges for each vertex
    std::vector<uint32_t> unique_counts(graph.num_nodes, 0);
    
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        
        if (start == end) continue; // No edges
        
        uint32_t last = UINT32_MAX; // Initialize to an impossible node ID
        uint32_t count = 0;
        
        for (uint32_t j = start; j < end; j++) {
            uint32_t neighbor = graph.col_idx[j];
            
            // Skip self-loops and duplicates
            if (neighbor != i && neighbor != last) {
                count++;
                last = neighbor;
            }
        }
        
        unique_counts[i] = count;
    }
    
    // Compute new row pointers based on unique counts
    std::vector<uint32_t> new_row_ptr(graph.num_nodes + 1);
    new_row_ptr[0] = 0;
    
    for (size_t i = 0; i < graph.num_nodes; i++) {
        new_row_ptr[i + 1] = new_row_ptr[i] + unique_counts[i];
    }
    
    // Allocate new column indices
    std::vector<uint32_t> new_col_idx(new_row_ptr.back());
    
    // Fill new column indices
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        uint32_t new_start = new_row_ptr[i];
        
        if (start == end) continue; // No edges
        
        uint32_t last = UINT32_MAX; // Initialize to an impossible node ID
        uint32_t pos = new_start;
        
        for (uint32_t j = start; j < end; j++) {
            uint32_t neighbor = graph.col_idx[j];
            
            // Skip self-loops and duplicates
            if (neighbor != i && neighbor != last) {
                new_col_idx[pos++] = neighbor;
                last = neighbor;
            }
        }
    }
    
    size_t removed_edges = (graph.col_idx.size() - new_col_idx.size()) / 2;
    
    // Update the graph
    graph.col_idx = std::move(new_col_idx);
    graph.row_ptr = std::move(new_row_ptr);
    
    if (verbose) {
        std::cout << "Removed " << removed_edges << " edges (self-loops and duplicates)" << std::endl;
    }
}

// Simple test function to validate the graph
void test_graph(const Graph& graph) {
    std::cout << "Testing graph integrity..." << std::endl;
    
    size_t total_edges = 0;
    for (size_t i = 0; i < graph.num_nodes; i++) {
        total_edges += graph.row_ptr[i + 1] - graph.row_ptr[i];
    }
    
    std::cout << "Total edges in CSR: " << total_edges << std::endl;
    
    // Check if all edges are valid
    bool valid = true;
    size_t checked_edges = 0;
    size_t total_to_check = std::min(size_t(1000), graph.num_nodes); // Limit check to 1000 nodes for performance
    
    for (size_t i = 0; i < total_to_check; i++) {
        for (size_t j = graph.row_ptr[i]; j < graph.row_ptr[i + 1]; j++) {
            uint32_t neighbor = graph.col_idx[j];
            checked_edges++;
            
            if (neighbor >= graph.num_nodes) {
                std::cout << "Invalid edge: " << i << " -> " << neighbor << std::endl;
                valid = false;
                break;
            }
            
            // Check if the reverse edge exists
            bool found = false;
            for (size_t k = graph.row_ptr[neighbor]; k < graph.row_ptr[neighbor + 1]; k++) {
                if (graph.col_idx[k] == i) {
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                std::cout << "Missing reverse edge: " << neighbor << " -> " << i << std::endl;
                valid = false;
                break;
            }
        }
        
        if (!valid) break;
    }
    
    if (valid) {
        std::cout << "Graph appears valid (checked " << checked_edges << " edges)" << std::endl;
    }
}

#endif // GRAPH_H
