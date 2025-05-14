#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <omp.h>
#include <utility>

// Double-Indexed (DI) representation for undirected graph
struct DIGraph {
    // Edge Index Array
    std::vector<uint32_t> src; // Source vertex array
    std::vector<uint32_t> dst; // Destination vertex array
    
    // Vertex Index Array
    std::vector<int32_t> str; // Starting position in edge arrays for each vertex
    std::vector<uint32_t> nei; // Number of edges for each vertex
    
    // Node ID mapping (if needed)
    std::unordered_map<uint64_t, uint32_t> node_map;
    std::vector<uint64_t> id_map;
    
    // Graph info
    size_t num_nodes = 0;
    size_t num_edges = 0; // This counts each undirected edge once
    
    // Initialize an empty graph with a specified number of nodes
    void init(size_t nodes) {
        num_nodes = nodes;
        // Initialize vertex index arrays
        str.resize(num_nodes, -1); // -1 indicates no edges
        nei.resize(num_nodes, 0);
        
        // Edge arrays will be populated as edges are added
        src.clear();
        dst.clear();
    }
    
    // Add an edge to the graph
    void add_edge(uint32_t from, uint32_t to) {
        // Check if nodes exist
        if (from >= num_nodes || to >= num_nodes) {
            return;
        }
        
        // Check if edge already exists (this is an O(nei[from]) operation)
        for (uint32_t i = 0; i < nei[from]; ++i) {
            if (str[from] != -1) {
                uint32_t edge_index = str[from] + i;
                if (edge_index < dst.size() && dst[edge_index] == to) {
                    return; // Edge already exists
                }
            }
        }
        
        // Add the edge (we'll sort and rebuild later)
        src.push_back(from);
        dst.push_back(to);
        
        // For undirected graph, add the reverse edge too
        src.push_back(to);
        dst.push_back(from);
        
        // Update edge count (only counting once for the undirected edge)
        num_edges++;
    }
    
    // Build the vertex index arrays after all edges have been added
    void build_index(int num_threads = 1, bool verbose = false) {
        if (verbose) {
            std::cout << "Building DI graph index..." << std::endl;
        }
        
        // First, create edge pairs for sorting
        std::vector<std::pair<std::pair<uint32_t, uint32_t>, uint32_t>> edges;
        edges.reserve(src.size());
        
        for (uint32_t i = 0; i < src.size(); i++) {
            edges.push_back({{src[i], dst[i]}, i});
        }
        
        // Sort edges by source vertex, then by destination vertex
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp single
            {
                std::sort(edges.begin(), edges.end());
            }
        }
        
        // Create new sorted arrays
        // Create new sorted arrays
        std::vector<uint32_t> new_src(src.size());
        std::vector<uint32_t> new_dst(dst.size());
        
        #pragma omp parallel for num_threads(num_threads)
        for (uint32_t i = 0; i < edges.size(); i++) {
            new_src[i] = edges[i].first.first;
            new_dst[i] = edges[i].first.second;
        }
        
        // Replace the old arrays with the sorted ones
        src = std::move(new_src);
        dst = std::move(new_dst);
        
        // Reset vertex index arrays
        #pragma omp parallel for num_threads(num_threads)
        for (size_t i = 0; i < num_nodes; i++) {
            str[i] = -1;
            nei[i] = 0;
        }
        
        // Count neighbors for each vertex (first pass)
        #pragma omp parallel for num_threads(num_threads)
        for (uint32_t i = 0; i < src.size(); i++) {
            #pragma omp atomic
            nei[src[i]]++;
        }
        
        // Compute starting positions (sequential, but fast)
        int32_t pos = 0;
        for (uint32_t v = 0; v < num_nodes; v++) {
            if (nei[v] > 0) {
            str[v] = pos;
            pos += nei[v];
            nei[v] = 0; // Reset for second pass
            }
        }
        
        // Build final neighbors count (second pass)
        #pragma omp parallel for num_threads(num_threads)
        for (uint32_t i = 0; i < src.size(); i++) {
            uint32_t v = src[i];
            #pragma omp atomic
            nei[v]++;
        }
        
        if (verbose) {
            std::cout << "DI graph index built with " << num_edges << " undirected edges" << std::endl;
        }
    }
};

// Remove self-loops and duplicate edges in parallel
void clean_graph_parallel(DIGraph& graph, int num_threads, bool verbose = false) {
    if (verbose) {
        std::cout << "Removing self-loops and duplicate edges..." << std::endl;
    }
    
    // We'll rebuild the graph from scratch
    std::vector<uint32_t> new_src;
    std::vector<uint32_t> new_dst;
    size_t new_edges = 0;
    
    // Process each vertex
    #pragma omp parallel num_threads(num_threads)
    {
        std::vector<uint32_t> local_src;
        std::vector<uint32_t> local_dst;
        size_t local_edges = 0;
        
        #pragma omp for schedule(dynamic, 1000)
        for (size_t v = 0; v < graph.num_nodes; v++) {
            if (graph.str[v] == -1) continue; // No edges for this vertex
            
            // Get all neighbors for this vertex
            std::vector<uint32_t> neighbors;
            for (uint32_t i = 0; i < graph.nei[v]; i++) {
                uint32_t edge_idx = graph.str[v] + i;
                uint32_t neighbor = graph.dst[edge_idx];
                
                // Skip self-loops
                if (neighbor != v) {
                    neighbors.push_back(neighbor);
                }
            }
            
            // Sort and remove duplicates
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
            
            // Add unique edges (but only if v < neighbor to avoid duplicating undirected edges)
            for (uint32_t neighbor : neighbors) {
                if (v < neighbor) {
                    local_src.push_back(v);
                    local_dst.push_back(neighbor);
                    local_src.push_back(neighbor);
                    local_dst.push_back(v);
                    local_edges++;
                }
            }
        }
        
        // Merge results
        #pragma omp critical
        {
            new_src.insert(new_src.end(), local_src.begin(), local_src.end());
            new_dst.insert(new_dst.end(), local_dst.begin(), local_dst.end());
            new_edges += local_edges;
        }
    }
    
    // Update the graph
    graph.src = std::move(new_src);
    graph.dst = std::move(new_dst);
    graph.num_edges = new_edges;
    
    // Rebuild the index
    graph.build_index(num_threads, verbose);
    
    if (verbose) {
        std::cout << "Clean graph has " << graph.num_edges << " undirected edges" << std::endl;
    }
}

// Simple test function to validate the graph
void test_graph(const DIGraph& graph) {
    std::cout << "Testing graph integrity..." << std::endl;
    
    size_t total_edges = 0;
    for (size_t i = 0; i < graph.num_nodes; i++) {
        total_edges += graph.nei[i];
    }
    
    std::cout << "Total edges in DI graph: " << total_edges << std::endl;
    
    // Check if all edges are valid
    bool valid = true;
    size_t checked_edges = 0;
    size_t total_to_check = std::min(size_t(1000), graph.num_nodes); // Limit check to 1000 nodes for performance
    
    for (size_t i = 0; i < total_to_check; i++) {
        if (graph.str[i] == -1) continue; // No edges for this vertex
        
        for (uint32_t j = 0; j < graph.nei[i]; j++) {
            uint32_t edge_idx = graph.str[i] + j;
            uint32_t neighbor = graph.dst[edge_idx];
            checked_edges++;
            
            if (neighbor >= graph.num_nodes) {
                std::cout << "Invalid edge: " << i << " -> " << neighbor << std::endl;
                valid = false;
                break;
            }
            
            // Check if the reverse edge exists
            bool found = false;
            if (graph.str[neighbor] != -1) {
                for (uint32_t k = 0; k < graph.nei[neighbor]; k++) {
                    uint32_t rev_edge_idx = graph.str[neighbor] + k;
                    if (graph.dst[rev_edge_idx] == i) {
                        found = true;
                        break;
                    }
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

// Utility function to print graph statistics
void print_graph_stats(const DIGraph& graph) {
    std::cout << "Graph Statistics:" << std::endl;
    std::cout << "Nodes: " << graph.num_nodes << std::endl;
    std::cout << "Edges: " << graph.num_edges << " (undirected)" << std::endl;
    
    // Calculate degree statistics
    if (graph.num_nodes > 0) {
        uint32_t max_degree = 0;
        uint32_t min_degree = UINT32_MAX;
        double avg_degree = 0.0;
        
        for (size_t i = 0; i < graph.num_nodes; i++) {
            uint32_t degree = graph.nei[i];
            max_degree = std::max(max_degree, degree);
            if (degree > 0) {
                min_degree = std::min(min_degree, degree);
            }
            avg_degree += degree;
        }
        
        avg_degree /= graph.num_nodes;
        
        std::cout << "Max degree: " << max_degree << std::endl;
        std::cout << "Min degree (non-zero): " << (min_degree == UINT32_MAX ? 0 : min_degree) << std::endl;
        std::cout << "Avg degree: " << avg_degree << std::endl;
    }
}

#endif // GRAPH_H
