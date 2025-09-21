#ifndef EDGE_EXTRACTOR_H
#define EDGE_EXTRACTOR_H

#include <vector>
#include <memory>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "../data_structures/graph.h"

// Utility class for extracting and processing edges from processed subgraphs
class EdgeExtractor {
public:
    // Extract all edges from completed subgraphs and map to original IDs
    static std::vector<std::pair<uint64_t, uint64_t>> extract_all_edges_from_subgraphs(
        const std::vector<std::shared_ptr<Graph>>& subgraphs) {
        
        std::vector<std::pair<uint64_t, uint64_t>> all_edges;
        
        for (const auto& subgraph : subgraphs) {
            auto subgraph_edges = extract_edges_from_single_subgraph(*subgraph);
            all_edges.insert(all_edges.end(), subgraph_edges.begin(), subgraph_edges.end());
        }
        
        // Remove duplicates (in case there are any overlaps)
        std::sort(all_edges.begin(), all_edges.end());
        all_edges.erase(std::unique(all_edges.begin(), all_edges.end()), all_edges.end());
        
        return all_edges;
    }
    
    // Extract edges from a single subgraph and convert to original IDs
    static std::vector<std::pair<uint64_t, uint64_t>> extract_edges_from_single_subgraph(
        const Graph& subgraph) {
        
        std::vector<std::pair<uint64_t, uint64_t>> edges;
        
        // Go through all edges in the subgraph
        for (uint32_t u = 0; u < subgraph.num_nodes; ++u) {
            for (uint32_t idx = subgraph.row_ptr[u]; idx < subgraph.row_ptr[u + 1]; ++idx) {
                uint32_t v = subgraph.col_idx[idx];
                
                // Only store each edge once (u < v)
                if (u < v) {
                    // Convert to original node IDs
                    uint64_t orig_u = subgraph.id_map[u];
                    uint64_t orig_v = subgraph.id_map[v];
                    
                    // Ensure consistent ordering
                    if (orig_u > orig_v) std::swap(orig_u, orig_v);
                    
                    edges.push_back({orig_u, orig_v});
                }
            }
        }
        
        return edges;
    }
    
    // Compare with original graph to find newly added edges
    static std::vector<std::pair<uint64_t, uint64_t>> find_newly_added_edges(
        const Graph& original_graph,
        const std::vector<std::shared_ptr<Graph>>& processed_subgraphs) {
        
        // Extract all edges from original graph
        std::unordered_set<uint64_t> original_edges_set;
        
        for (uint32_t u = 0; u < original_graph.num_nodes; ++u) {
            for (uint32_t idx = original_graph.row_ptr[u]; idx < original_graph.row_ptr[u + 1]; ++idx) {
                uint32_t v = original_graph.col_idx[idx];
                
                if (u < v) {
                    // Convert to original IDs
                    uint64_t orig_u = original_graph.id_map[u];
                    uint64_t orig_v = original_graph.id_map[v];
                    if (orig_u > orig_v) std::swap(orig_u, orig_v);
                    
                    // Pack into single uint64_t for fast lookup
                    uint64_t edge_key = (orig_u << 32) | orig_v;
                    original_edges_set.insert(edge_key);
                }
            }
        }
        
        // Extract all edges from processed subgraphs
        auto all_processed_edges = extract_all_edges_from_subgraphs(processed_subgraphs);
        
        // Find newly added edges
        std::vector<std::pair<uint64_t, uint64_t>> newly_added_edges;
        
        for (const auto& [u, v] : all_processed_edges) {
            uint64_t edge_key = (u << 32) | v;
            
            if (original_edges_set.find(edge_key) == original_edges_set.end()) {
                newly_added_edges.push_back({u, v});
            }
        }
        
        return newly_added_edges;
    }
    
    // Write edges to TSV file
    static void write_edges_to_tsv(const std::vector<std::pair<uint64_t, uint64_t>>& edges,
                                   const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
            return;
        }
        
        // Write header
        file << "node1\tnode2\n";
        
        // Write edges
        for (const auto& [u, v] : edges) {
            file << u << "\t" << v << "\n";
        }
        
        std::cout << "Wrote " << edges.size() << " edges to " << filename << std::endl;
    }
    
    // Write edges to TSV file without header
    static void write_edges_to_tsv_no_header(const std::vector<std::pair<uint64_t, uint64_t>>& edges,
                                            const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
            return;
        }
        
        // Write edges without header
        for (const auto& [u, v] : edges) {
            file << u << "\t" << v << "\n";
        }
        
        std::cout << "Wrote " << edges.size() << " edges to " << filename << std::endl;
    }

    // Fetch the compressed ids of the newly added edges
    static std::vector<std::pair<uint32_t, uint32_t>> get_compressed_newly_added_edges(
        const Graph& original_graph,
        const std::vector<std::pair<uint64_t, uint64_t>>& newly_added_edges_raw,
        bool verbose = false) {
        
        // Convert original node IDs back to internal indices
        std::vector<std::pair<uint32_t, uint32_t>> newly_added_edges;
        newly_added_edges.reserve(newly_added_edges_raw.size());

        for (const auto& edge : newly_added_edges_raw) {
            // Convert from original IDs to internal indices, adding nodes if missing
            auto it_u = original_graph.node_map.find(edge.first);
            auto it_v = original_graph.node_map.find(edge.second);
            
            uint32_t internal_u, internal_v;
            
            // Handle node u
            if (it_u != original_graph.node_map.end()) {
                internal_u = it_u->second;
            } else {
                if (verbose) {
                    std::cout << "Adding missing node " << edge.first << " to graph" << std::endl;
                }
                // Note: This assumes original_graph has an add_node method
                // If not, this may need to be handled differently
                internal_u = const_cast<Graph&>(original_graph).add_node(edge.first);
            }
            
            // Handle node v
            if (it_v != original_graph.node_map.end()) {
                internal_v = it_v->second;
            } else {
                if (verbose) {
                    std::cout << "Adding missing node " << edge.second << " to graph" << std::endl;
                }
                internal_v = const_cast<Graph&>(original_graph).add_node(edge.second);
            }
            
            newly_added_edges.emplace_back(internal_u, internal_v);
        }
        
        return newly_added_edges;
    }
};

#endif // EDGE_EXTRACTOR_H
