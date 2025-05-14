#ifndef GRAPH_MERGER_H
#define GRAPH_MERGER_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "../data_structures/graph.h"

// Helper struct for hashing std::pair
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1);
    }
};


class GraphMerger {
public:
    // Append a DIGraph to an existing edge list file
    static bool append_graph_to_file(
        const DIGraph& graph,
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
        
        // Use a set to track edges we've already processed to avoid duplicates
        std::unordered_set<std::pair<uint32_t, uint32_t>, 
            pair_hash> processed_edges;
        
        // For each edge in the graph
        for (size_t i = 0; i < graph.src.size(); i++) {
            uint32_t node_id = graph.src[i];
            uint32_t neighbor_id = graph.dst[i];
            
            // Create an ordered pair to ensure we only process each edge once
            uint32_t src = std::min(node_id, neighbor_id);
            uint32_t dst = std::max(node_id, neighbor_id);
            
            // Check if we've already processed this edge
            auto edge_pair = std::make_pair(src, dst);
            if (processed_edges.find(edge_pair) != processed_edges.end()) {
                continue;
            }
            
            // Mark as processed
            processed_edges.insert(edge_pair);
            
            // Map to original IDs
            uint64_t original_src_id = graph.id_map[src];
            uint64_t original_dst_id = graph.id_map[dst];
            
            // Write edge to file
            out << original_src_id << "\t" << original_dst_id << "\n";
            edges_appended++;
        }
        
        out.close();
        
        if (verbose) {
            std::cout << "Successfully appended " << edges_appended << " edges to the output file" << std::endl;
        }
        
        return true;
    }
    
    // Function that takes a modified clustered graph and appends it to an unclustered graph file
    static bool merge_clustered_to_unclustered(
        const DIGraph& clustered_graph,
        const std::string& unclustered_file,
        const std::string& output_file,
        bool verbose = false) {
        
        return append_graph_to_file(clustered_graph, unclustered_file, output_file, verbose);
    }
};
#endif // GRAPH_MERGER_H
