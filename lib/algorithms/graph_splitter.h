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
    static std::tuple<CSRGraph, CSRGraph, std::unordered_map<uint32_t, uint32_t>, std::unordered_map<uint32_t, uint32_t>> 
    split_by_clustering(
        const CSRGraph& original_graph, 
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
    static std::tuple<CSRGraph, CSRGraph, std::unordered_map<uint32_t, uint32_t>, std::unordered_map<uint32_t, uint32_t>> 
    split_graph_single_pass(
        const CSRGraph& original_graph,
        const std::unordered_set<uint32_t>& clustered_nodes,
        bool verbose) {
        
        // Create the clustered subgraph
        CSRGraph clustered_graph;
        
        // Create mapping from original node IDs to new node IDs for clustered subgraph
        std::unordered_map<uint32_t, uint32_t> clustered_node_map;
        
        uint32_t clustered_idx = 0;
        for (uint32_t node_id : clustered_nodes) {
            clustered_node_map[node_id] = clustered_idx++;
        }
        
        // Track nodes that have edges in the complement graph
        std::unordered_set<uint32_t> complement_nodes_with_edges;
        
        // First pass: count edges for both subgraphs
        std::vector<uint32_t> clustered_degree(clustered_nodes.size(), 0);
        std::vector<uint32_t> complement_degree_temp;
        complement_degree_temp.resize(original_graph.num_nodes, 0);
        
        size_t clustered_edge_count = 0;
        size_t complement_edge_count = 0;
        
        for (uint32_t node_id = 0; node_id < original_graph.num_nodes; ++node_id) {
            // Get original node edges
            for (uint32_t i = original_graph.row_ptr[node_id]; 
                 i < original_graph.row_ptr[node_id + 1]; ++i) {
                uint32_t neighbor = original_graph.col_idx[i];
                
                // Check if both endpoints are in clustered nodes
                bool node_in_cluster = (clustered_nodes.find(node_id) != clustered_nodes.end());
                bool neighbor_in_cluster = (clustered_nodes.find(neighbor) != clustered_nodes.end());
                
                if (node_in_cluster && neighbor_in_cluster) {
                    // This edge belongs to the clustered subgraph
                    uint32_t new_node_id = clustered_node_map[node_id];
                    clustered_degree[new_node_id]++;
                    clustered_edge_count++;
                } else {
                    // This edge belongs to the complement subgraph
                    complement_degree_temp[node_id]++;
                    complement_edge_count++;
                    
                    // Track nodes that will have edges in the complement graph
                    complement_nodes_with_edges.insert(node_id);
                    complement_nodes_with_edges.insert(neighbor);
                }
            }
        }
        
        // Set final edge counts (divide by 2 because we counted each edge twice in undirected graph)
        clustered_graph.num_edges = clustered_edge_count / 2;
        
        if (verbose) {
            std::cout << "Complement nodes with edges: " << complement_nodes_with_edges.size() << std::endl;
        }
        
        // Create the complement graph
        CSRGraph complement_graph;
        
        // Create mapping for complement graph nodes
        std::unordered_map<uint32_t, uint32_t> complement_node_map;
        std::vector<uint32_t> complement_reverse_map(complement_nodes_with_edges.size());
        
        uint32_t complement_idx = 0;
        for (uint32_t node_id : complement_nodes_with_edges) {
            complement_node_map[node_id] = complement_idx;
            complement_reverse_map[complement_idx] = node_id;
            complement_idx++;
        }
        
        // Prepare the clustered graph
        clustered_graph.num_nodes = clustered_nodes.size();
        clustered_graph.row_ptr.resize(clustered_graph.num_nodes + 1, 0);
        
        // Compute row pointers for clustered graph
        for (uint32_t i = 0; i < clustered_graph.num_nodes; ++i) {
            clustered_graph.row_ptr[i + 1] = clustered_graph.row_ptr[i] + clustered_degree[i];
        }
        
        // Prepare the complement graph
        complement_graph.num_nodes = complement_nodes_with_edges.size();
        complement_graph.num_edges = complement_edge_count / 2;
        complement_graph.row_ptr.resize(complement_graph.num_nodes + 1, 0);
        
        // Create properly sized degree array for complement graph
        std::vector<uint32_t> complement_degree(complement_graph.num_nodes, 0);
        
        // Update the complement degree array with the proper indices
        for (uint32_t new_idx = 0; new_idx < complement_graph.num_nodes; ++new_idx) {
            uint32_t orig_node_id = complement_reverse_map[new_idx];
            complement_degree[new_idx] = complement_degree_temp[orig_node_id];
        }
        
        // Compute row pointers for complement graph
        for (uint32_t i = 0; i < complement_graph.num_nodes; ++i) {
            complement_graph.row_ptr[i + 1] = complement_graph.row_ptr[i] + complement_degree[i];
        }
        
        // Reset degree arrays for filling column indices
        std::fill(clustered_degree.begin(), clustered_degree.end(), 0);
        std::fill(complement_degree.begin(), complement_degree.end(), 0);
        
        // Allocate column indices
        clustered_graph.col_idx.resize(clustered_graph.row_ptr.back());
        complement_graph.col_idx.resize(complement_graph.row_ptr.back());
        
        // Second pass: fill column indices for both subgraphs
        for (uint32_t node_id = 0; node_id < original_graph.num_nodes; ++node_id) {
            // Get original node edges
            for (uint32_t i = original_graph.row_ptr[node_id]; 
                 i < original_graph.row_ptr[node_id + 1]; ++i) {
                uint32_t neighbor = original_graph.col_idx[i];
                
                // Check if both endpoints are in clustered nodes
                bool node_in_cluster = (clustered_nodes.find(node_id) != clustered_nodes.end());
                bool neighbor_in_cluster = (clustered_nodes.find(neighbor) != clustered_nodes.end());
                
                if (node_in_cluster && neighbor_in_cluster) {
                    // This edge belongs to the clustered subgraph
                    uint32_t new_node_id = clustered_node_map[node_id];
                    uint32_t new_neighbor_id = clustered_node_map[neighbor];
                    
                    uint32_t pos = clustered_graph.row_ptr[new_node_id] + clustered_degree[new_node_id];
                    clustered_graph.col_idx[pos] = new_neighbor_id;
                    clustered_degree[new_node_id]++;
                } else {
                    // This edge belongs to the complement subgraph
                    auto it_node = complement_node_map.find(node_id);
                    auto it_neighbor = complement_node_map.find(neighbor);
                    
                    // Both nodes should be in the map since we added them earlier
                    if (it_node != complement_node_map.end() && it_neighbor != complement_node_map.end()) {
                        uint32_t new_node_id = it_node->second;
                        uint32_t new_neighbor_id = it_neighbor->second;
                        
                        uint32_t pos = complement_graph.row_ptr[new_node_id] + complement_degree[new_node_id];
                        complement_graph.col_idx[pos] = new_neighbor_id;
                        complement_degree[new_node_id]++;
                    }
                }
            }
        }
        
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
