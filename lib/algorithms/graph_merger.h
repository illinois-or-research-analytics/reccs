#ifndef GRAPH_MERGER_H
#define GRAPH_MERGER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "../data_structures/graph.h"

class GraphMerger {
public:
    // Append a CSRGraph to an existing edge list file
    static bool append_graph_to_file(
        const CSRGraph& graph,
        const std::string& existing_file,
        const std::string& output_file,
        bool verbose = false) {
        
        if (verbose) {
            std::cout << "Appending graph to existing edge list..." << std::endl;
            std::cout << "Graph to append: " << graph.num_nodes << " nodes, " 
                      << graph.num_edges << " edges" << std::endl;
            std::cout << "Existing file: " << existing_file << std::endl;
            std::cout << "Output file: " << output_file << std::endl;
        }
        
        // First, copy the existing file to the output file
        std::ifstream src(existing_file, std::ios::binary);
        if (!src.is_open()) {
            std::cerr << "Failed to open existing file: " << existing_file << std::endl;
            return false;
        }
        
        std::ofstream dst(output_file, std::ios::binary);
        if (!dst.is_open()) {
            std::cerr << "Failed to open output file: " << output_file << std::endl;
            src.close();
            return false;
        }
        
        dst << src.rdbuf();
        src.close();
        
        // Now append the graph edges to the output file
        // The file is already opened in binary mode, so we'll close and reopen in text mode
        dst.close();
        std::ofstream out(output_file, std::ios::app);  // Open in append mode
        
        if (!out.is_open()) {
            std::cerr << "Failed to reopen output file for appending: " << output_file << std::endl;
            return false;
        }
        
        // Count the edges we're appending
        size_t edges_appended = 0;
        
        // For each node in the graph
        for (uint32_t node_id = 0; node_id < graph.num_nodes; ++node_id) {
            uint64_t original_node_id = graph.id_map[node_id];
            
            // For each neighbor of this node
            for (uint32_t i = graph.row_ptr[node_id]; i < graph.row_ptr[node_id + 1]; ++i) {
                uint32_t neighbor_id = graph.col_idx[i];
                
                // Only process each edge once (where source < target)
                if (node_id < neighbor_id) {
                    uint64_t original_neighbor_id = graph.id_map[neighbor_id];
                    
                    // Write edge to file
                    out << original_node_id << "\t" << original_neighbor_id << "\n";
                    edges_appended++;
                }
            }
        }
        
        out.close();
        
        if (verbose) {
            std::cout << "Successfully appended " << edges_appended << " edges to the output file" << std::endl;
        }
        
        return true;
    }
    
    // Function that takes a modified clustered graph and appends it to an unclustered graph file
    static bool merge_clustered_to_unclustered(
        const CSRGraph& clustered_graph,
        const std::string& unclustered_file,
        const std::string& output_file,
        bool verbose = false) {
        
        return append_graph_to_file(clustered_graph, unclustered_file, output_file, verbose);
    }
};

#endif // GRAPH_MERGER_H
