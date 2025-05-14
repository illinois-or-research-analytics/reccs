#ifndef GRAPH_SPLITTER_H
#define GRAPH_SPLITTER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <tuple>
#include "../data_structures/graph.h"
#include "../data_structures/clustering.h"

class GraphSplitter {
public:
    // Split a graph based on clustering into:
    // 1. Nodes that are part of non-singleton clusters and their internal edges
    // 2. The complement graph containing all edges not in the clustered subgraph (excluding isolated nodes)
    // Returns: (clustered_graph, complement_graph, clustered_node_map, complement_node_map)
    static std::tuple<DIGraph, DIGraph, std::unordered_map<uint32_t, uint32_t>, std::unordered_map<uint32_t, uint32_t>> 
    split_by_clustering(
        const DIGraph& original_graph, 
        const Clustering& clustering,
        bool verbose = false) {
        
        if (verbose) {
            std::cout << "Splitting graph based on clustering..." << std::endl;
        }
        
        // First, identify singleton clusters
        std::unordered_set<uint32_t> singleton_cluster_indices;
        for (uint32_t cluster_idx = 0; cluster_idx < clustering.cluster_nodes.size(); ++cluster_idx) {
            if (clustering.cluster_nodes[cluster_idx].size() == 1) {
                singleton_cluster_indices.insert(cluster_idx);
            }
        }
        
        if (verbose) {
            std::cout << "Found " << singleton_cluster_indices.size() 
                      << " singleton clusters out of " 
                      << clustering.get_non_empty_cluster_count() << " total clusters" << std::endl;
        }
        
        // Identify nodes in non-singleton clusters
        std::unordered_set<uint32_t> clustered_nodes;
        
        for (uint32_t node_id = 0; node_id < original_graph.num_nodes; ++node_id) {
            if (node_id < clustering.node_to_cluster_idx.size() && 
                clustering.node_to_cluster_idx[node_id] != UINT32_MAX) {
                
                uint32_t cluster_idx = clustering.node_to_cluster_idx[node_id];
                if (singleton_cluster_indices.find(cluster_idx) == singleton_cluster_indices.end()) {
                    // Node is in a non-singleton cluster
                    clustered_nodes.insert(node_id);
                }
            }
        }
        
        if (verbose) {
            std::cout << "Nodes in non-singleton clusters: " << clustered_nodes.size() << std::endl;
            std::cout << "Nodes not in non-singleton clusters: " 
                      << (original_graph.num_nodes - clustered_nodes.size()) << std::endl;
        }
        
        // Single-pass graph splitting
        return split_graph_single_pass(original_graph, clustered_nodes, verbose);
    }

private:
    // Split graph in a single pass and return both subgraphs and their node mappings
    static std::tuple<DIGraph, DIGraph, std::unordered_map<uint32_t, uint32_t>, std::unordered_map<uint32_t, uint32_t>> 
    split_graph_single_pass(
        const DIGraph& original_graph,
        const std::unordered_set<uint32_t>& clustered_nodes,
        bool verbose) {
        
        // Create the clustered subgraph
        DIGraph clustered_graph;
        
        // Create mapping from original node IDs to new node IDs for clustered subgraph
        std::unordered_map<uint32_t, uint32_t> clustered_node_map;
        
        uint32_t clustered_idx = 0;
        for (uint32_t node_id : clustered_nodes) {
            clustered_node_map[node_id] = clustered_idx++;
        }
        
        // Track nodes that have edges in the complement graph
        std::unordered_set<uint32_t> complement_nodes_with_edges;
        
        // Initialize the subgraphs
        clustered_graph.init(clustered_nodes.size());
        clustered_graph.num_nodes = clustered_nodes.size();
        clustered_graph.num_edges = 0;
        
        // First pass: identify edges for both subgraphs
        std::vector<std::pair<uint32_t, uint32_t>> clustered_edges;
        std::vector<std::pair<uint32_t, uint32_t>> complement_edges;
        
        // Process each edge in the original graph
        for (uint32_t i = 0; i < original_graph.src.size(); i++) {
            uint32_t node_id = original_graph.src[i];
            uint32_t neighbor = original_graph.dst[i];
            
            // Check if both endpoints are in clustered nodes
            bool node_in_cluster = (clustered_nodes.find(node_id) != clustered_nodes.end());
            bool neighbor_in_cluster = (clustered_nodes.find(neighbor) != clustered_nodes.end());
            
            if (node_in_cluster && neighbor_in_cluster) {
                // This edge belongs to the clustered subgraph
                if (node_id < neighbor) { // Count each undirected edge only once
                    uint32_t new_node_id = clustered_node_map[node_id];
                    uint32_t new_neighbor_id = clustered_node_map[neighbor];
                    clustered_edges.emplace_back(new_node_id, new_neighbor_id);
                }
            } else {
                // This edge belongs to the complement subgraph
                if (node_id < neighbor) { // Count each undirected edge only once
                    complement_edges.emplace_back(node_id, neighbor);
                    
                    // Track nodes that will have edges in the complement graph
                    complement_nodes_with_edges.insert(node_id);
                    complement_nodes_with_edges.insert(neighbor);
                }
            }
        }
        
        // Create mapping for complement graph nodes
        std::unordered_map<uint32_t, uint32_t> complement_node_map;
        std::vector<uint32_t> complement_reverse_map(complement_nodes_with_edges.size());
        
        uint32_t complement_idx = 0;
        for (uint32_t node_id : complement_nodes_with_edges) {
            complement_node_map[node_id] = complement_idx;
            complement_reverse_map[complement_idx] = node_id;
            complement_idx++;
        }
        
        // Create and initialize the complement graph
        DIGraph complement_graph;
        complement_graph.init(complement_nodes_with_edges.size());
        complement_graph.num_nodes = complement_nodes_with_edges.size();
        complement_graph.num_edges = 0;
        
        // Add edges to the clustered graph
        for (const auto& [src, dst] : clustered_edges) {
            // Add in both directions (for undirected graph)
            clustered_graph.src.push_back(src);
            clustered_graph.dst.push_back(dst);
            clustered_graph.src.push_back(dst);
            clustered_graph.dst.push_back(src);
            clustered_graph.num_edges++; // Count as one undirected edge
        }
        
        // Add edges to the complement graph
        for (const auto& [orig_src, orig_dst] : complement_edges) {
            uint32_t new_src = complement_node_map[orig_src];
            uint32_t new_dst = complement_node_map[orig_dst];
            
            // Add in both directions (for undirected graph)
            complement_graph.src.push_back(new_src);
            complement_graph.dst.push_back(new_dst);
            complement_graph.src.push_back(new_dst);
            complement_graph.dst.push_back(new_src);
            complement_graph.num_edges++; // Count as one undirected edge
        }
        
        // Build indices for both graphs
        clustered_graph.build_index();
        complement_graph.build_index();
        
        // Build the ID maps for clustered graph
        clustered_graph.id_map.resize(clustered_graph.num_nodes);
        for (const auto& [orig_node_id, new_node_id] : clustered_node_map) {
            clustered_graph.id_map[new_node_id] = original_graph.id_map[orig_node_id];
            clustered_graph.node_map[original_graph.id_map[orig_node_id]] = new_node_id;
        }
        
        // Build the ID maps for complement graph
        complement_graph.id_map.resize(complement_graph.num_nodes);
        for (uint32_t new_node_id = 0; new_node_id < complement_graph.num_nodes; ++new_node_id) {
            uint32_t orig_node_id = complement_reverse_map[new_node_id];
            complement_graph.id_map[new_node_id] = original_graph.id_map[orig_node_id];
            complement_graph.node_map[original_graph.id_map[orig_node_id]] = new_node_id;
        }
        
        if (verbose) {
            std::cout << "Split graph in a single pass" << std::endl;
            std::cout << "Clustered subgraph: " << clustered_graph.num_nodes << " nodes, " 
                      << clustered_graph.num_edges << " edges" << std::endl;
            std::cout << "Complement graph: " << complement_graph.num_nodes << " nodes, " 
                      << complement_graph.num_edges << " edges" << std::endl;
        }
        
        return {clustered_graph, complement_graph, clustered_node_map, complement_node_map};
    }
};

#endif // GRAPH_SPLITTER_H
