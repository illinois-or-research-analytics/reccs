#ifndef DEGREE_ENFORCER_H
#define DEGREE_ENFORCER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <iostream>
#include "../data_structures/graph.h"
#include "../data_structures/clustering.h"

class DegreeEnforcer {
public:
    // Enforce minimum cluster degrees based on original graph
    static CSRGraph enforce_min_cluster_degrees(
        const CSRGraph& sbm_graph,                  // Generated SBM graph 
        const CSRGraph& original_graph,             // Original graph for reference
        const Clustering& clustering,               // Clustering information
        const std::unordered_map<uint32_t, uint32_t>& original_to_sbm_node_map, // Mapping
        uint32_t minimum_degree_floor = 1,          // Minimum degree floor (override 0)
        bool verbose = false) {
        
        if (verbose) {
            std::cout << "Enforcing minimum cluster degrees..." << std::endl;
        }
        
        // Create a new graph by collecting edges
        std::vector<std::pair<uint64_t, uint64_t>> all_edges;
        std::unordered_map<uint64_t, uint32_t> all_nodes;  // Maps original ID to new index
        
        // Collect existing edges from SBM graph
        for (uint32_t node_id = 0; node_id < sbm_graph.num_nodes; ++node_id) {
            uint64_t src_id = sbm_graph.id_map[node_id];
            all_nodes[src_id] = 0;  // We'll assign indices later
            
            for (uint32_t i = sbm_graph.row_ptr[node_id]; i < sbm_graph.row_ptr[node_id + 1]; ++i) {
                uint32_t neighbor = sbm_graph.col_idx[i];
                if (node_id < neighbor) {  // Only count each edge once
                    uint64_t dst_id = sbm_graph.id_map[neighbor];
                    all_nodes[dst_id] = 0;  // We'll assign indices later
                    all_edges.push_back(std::make_pair(src_id, dst_id));
                }
            }
        }
        
        if (verbose) {
            std::cout << "Collected " << all_edges.size() << " edges from SBM graph" << std::endl;
        }
        
        // 1. Compute minimum degree per cluster in the original graph
        std::unordered_map<uint32_t, uint32_t> cluster_min_degrees;
        
        for (uint32_t cluster_idx = 0; cluster_idx < clustering.cluster_nodes.size(); ++cluster_idx) {
            // Skip empty clusters
            if (clustering.cluster_nodes[cluster_idx].empty()) {
                continue;
            }
            
            // Get the nodes in this cluster
            const auto& nodes = clustering.cluster_nodes[cluster_idx];
            
            // Find the minimum degree within the cluster
            uint32_t min_degree = UINT32_MAX;
            
            for (uint32_t node_id : nodes) {
                if (node_id >= original_graph.num_nodes) {
                    continue;  // Skip invalid nodes
                }
                
                // Count edges to nodes in the same cluster
                uint32_t internal_degree = 0;
                
                for (uint32_t i = original_graph.row_ptr[node_id]; 
                     i < original_graph.row_ptr[node_id + 1]; ++i) {
                    uint32_t neighbor = original_graph.col_idx[i];
                    
                    // Check if neighbor is in the same cluster
                    if (nodes.find(neighbor) != nodes.end()) {
                        internal_degree++;
                    }
                }
                
                min_degree = std::min(min_degree, internal_degree);
            }
            
            // Store the minimum degree for this cluster (with a floor)
            if (min_degree != UINT32_MAX) {
                // Apply the floor - never go below minimum_degree_floor
                min_degree = std::max(min_degree, minimum_degree_floor);
                cluster_min_degrees[cluster_idx] = min_degree;
                
                if (verbose) {
                    std::cout << "  Cluster " << cluster_idx << ": minimum degree = " 
                              << min_degree << std::endl;
                }
            }
        }
        
        // 2. Map nodes to clusters
        // Create a mapping from original node IDs to cluster IDs
        std::unordered_map<uint64_t, uint32_t> node_id_to_cluster;
        std::unordered_map<uint32_t, std::vector<uint64_t>> clusters_to_nodes;
        
        for (const auto& [orig_node, sbm_node] : original_to_sbm_node_map) {
            if (orig_node >= clustering.node_to_cluster_idx.size()) {
                continue;  // Skip invalid nodes
            }
            
            uint32_t cluster_idx = clustering.node_to_cluster_idx[orig_node];
            if (cluster_idx != UINT32_MAX) {
                // Get the original node ID from the SBM graph
                if (sbm_node < sbm_graph.num_nodes) {
                    uint64_t orig_id = sbm_graph.id_map[sbm_node];
                    node_id_to_cluster[orig_id] = cluster_idx;
                    clusters_to_nodes[cluster_idx].push_back(orig_id);
                }
            }
        }
        
        if (verbose) {
            std::cout << "  Mapped " << node_id_to_cluster.size() << " nodes to clusters" << std::endl;
            std::cout << "  Found " << clusters_to_nodes.size() << " clusters in SBM graph" << std::endl;
        }
        
        // 3. Compute current internal degrees and build adjacency structures
        std::unordered_map<uint64_t, uint32_t> node_internal_degrees;
        std::unordered_map<uint64_t, std::unordered_set<uint64_t>> adjacency_sets;

        // Initialize structures
        for (const auto& [cluster_idx, nodes] : clusters_to_nodes) {
            for (uint64_t node_id : nodes) {
                node_internal_degrees[node_id] = 0;
                adjacency_sets[node_id] = std::unordered_set<uint64_t>();
            }
        }

        // Build adjacency sets and count internal degrees
        for (const auto& edge : all_edges) {
            uint64_t src_id = edge.first;
            uint64_t dst_id = edge.second;
            
            // Add to adjacency sets
            adjacency_sets[src_id].insert(dst_id);
            adjacency_sets[dst_id].insert(src_id);
            
            // Count internal degrees
            auto src_cluster_it = node_id_to_cluster.find(src_id);
            auto dst_cluster_it = node_id_to_cluster.find(dst_id);
            
            if (src_cluster_it != node_id_to_cluster.end() && 
                dst_cluster_it != node_id_to_cluster.end() && 
                src_cluster_it->second == dst_cluster_it->second) {
                
                // This is an internal edge
                node_internal_degrees[src_id]++;
                node_internal_degrees[dst_id]++;
            }
        }

        // 4. Add edges to ensure minimum degree requirements
        // Setup random number generator
        std::random_device rd;
        std::mt19937 gen(rd());

        // Track added edges for reporting
        std::atomic<size_t> edges_added(0);

        // Convert map to vector for parallelization
        std::vector<std::pair<uint32_t, std::vector<uint64_t>>> clusters_vector;
        for (const auto& cluster_pair : clusters_to_nodes) {
            clusters_vector.push_back(cluster_pair);
        }

        // Process each cluster in parallel
        #pragma omp parallel
        {
            // Create a thread-local random generator
            std::mt19937 local_gen(rd() + omp_get_thread_num());
            
            // Process clusters in parallel
            #pragma omp for
            for (size_t i = 0; i < clusters_vector.size(); i++) {
                uint32_t cluster_idx = clusters_vector[i].first;
                const std::vector<uint64_t>& nodes = clusters_vector[i].second;
                
                // Skip clusters with no minimum degree requirement
                if (cluster_min_degrees.find(cluster_idx) == cluster_min_degrees.end()) {
                    continue;
                }
                
                uint32_t min_degree = cluster_min_degrees[cluster_idx];
                std::vector<std::pair<uint64_t, uint64_t>> local_edges_to_add;
                
                // For each node, ensure it has at least min_degree cluster neighbors
                for (uint64_t node_id : nodes) {
                    // Thread-local copies to avoid race conditions
                    auto local_internal_degree = node_internal_degrees[node_id];
                    auto local_adjacency = adjacency_sets[node_id];
                    
                    // Skip nodes that already have sufficient degree
                    if (local_internal_degree >= min_degree) {
                        continue;
                    }
                    
                    // How many more edges do we need?
                    uint32_t edges_needed = min_degree - local_internal_degree;
                    
                    // Find potential neighbors efficiently
                    std::vector<uint64_t> potential_neighbors;
                    for (uint64_t other_node : nodes) {
                        // Skip self-loops and existing neighbors
                        if (other_node == node_id || local_adjacency.count(other_node) > 0) {
                            continue;
                        }
                        
                        potential_neighbors.push_back(other_node);
                    }
                    
                    // If we don't have enough potential neighbors, we can't satisfy the requirement
                    if (potential_neighbors.size() < edges_needed) {
                        #pragma omp critical
                        {
                            if (verbose) {
                                std::cout << "  Warning: Node " << node_id << " in cluster " << cluster_idx
                                        << " needs " << edges_needed << " more edges, but only "
                                        << potential_neighbors.size() << " potential neighbors available" << std::endl;
                            }
                        }
                        
                        // Add as many edges as we can
                        edges_needed = potential_neighbors.size();
                    }
                    
                    // Randomly select neighbors to add
                    std::shuffle(potential_neighbors.begin(), potential_neighbors.end(), local_gen);
                    
                    for (uint32_t j = 0; j < edges_needed && j < potential_neighbors.size(); ++j) {
                        uint64_t new_neighbor = potential_neighbors[j];
                        
                        // Add edge to our local collection
                        local_edges_to_add.push_back(std::make_pair(
                            std::min(node_id, new_neighbor), 
                            std::max(node_id, new_neighbor)));
                        
                        // Update local adjacency sets
                        local_adjacency.insert(new_neighbor);
                    }
                    
                    // Store updated adjacency for later merging
                    #pragma omp critical
                    {
                        adjacency_sets[node_id].insert(local_adjacency.begin(), local_adjacency.end());
                    }
                }
                
                // Add local edges to global collection
                #pragma omp critical
                {
                    all_edges.insert(all_edges.end(), local_edges_to_add.begin(), local_edges_to_add.end());
                    edges_added += local_edges_to_add.size();
                }
            }
        }
        
        if (verbose) {
            std::cout << "Added " << edges_added << " edges to enforce minimum cluster degrees" << std::endl;
        }
        
        // 5. Build new CSRGraph from collected edges
        // First, assign indices to all nodes
        uint32_t next_idx = 0;
        for (auto& entry : all_nodes) {
            entry.second = next_idx++;
        }
        
        CSRGraph new_graph;
        new_graph.num_nodes = all_nodes.size();
        
        // Create ID mappings
        new_graph.id_map.resize(new_graph.num_nodes);
        for (const auto& [orig_id, idx] : all_nodes) {
            new_graph.id_map[idx] = orig_id;
            new_graph.node_map[orig_id] = idx;
        }
        
        // Count edges per node to create row pointers
        std::vector<uint32_t> node_edge_counts(new_graph.num_nodes, 0);
        
        for (const auto& edge : all_edges) {
            uint32_t src_idx = all_nodes[edge.first];
            uint32_t dst_idx = all_nodes[edge.second];
            
            node_edge_counts[src_idx]++;
            node_edge_counts[dst_idx]++;
        }
        
        // Create row pointers
        new_graph.row_ptr.resize(new_graph.num_nodes + 1, 0);
        for (uint32_t i = 0; i < new_graph.num_nodes; ++i) {
            new_graph.row_ptr[i + 1] = new_graph.row_ptr[i] + node_edge_counts[i];
        }
        
        // Reset edge counts for filling column indices
        std::fill(node_edge_counts.begin(), node_edge_counts.end(), 0);
        
        // Create column indices array
        new_graph.col_idx.resize(new_graph.row_ptr.back());
        
        // Fill column indices
        for (const auto& edge : all_edges) {
            uint32_t src_idx = all_nodes[edge.first];
            uint32_t dst_idx = all_nodes[edge.second];
            
            // Add edge in both directions
            uint32_t pos1 = new_graph.row_ptr[src_idx] + node_edge_counts[src_idx];
            new_graph.col_idx[pos1] = dst_idx;
            node_edge_counts[src_idx]++;
            
            uint32_t pos2 = new_graph.row_ptr[dst_idx] + node_edge_counts[dst_idx];
            new_graph.col_idx[pos2] = src_idx;
            node_edge_counts[dst_idx]++;
        }
        
        // Set edge count
        new_graph.num_edges = all_edges.size();
        
        if (verbose) {
            std::cout << "Created new graph with " << new_graph.num_nodes << " nodes and " 
                      << new_graph.num_edges << " edges" << std::endl;
        }
        
        return new_graph;
    }
};

#endif // DEGREE_ENFORCER_H
