#ifndef CLUSTER_IO_H
#define CLUSTER_IO_H

#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <chrono>
#include <iomanip>
#include "../data_structures/clustering.h"
#include "../data_structures/graph.h"
#include "mapped_file.h"

// Helper struct for parallel loading
struct ClusteringParseResult {
    std::unordered_map<std::string, std::unordered_set<uint64_t>> cluster_original_nodes;
};

// Parse a range of the file in parallel
ClusteringParseResult parse_clustering_chunk(const char* data, size_t begin, size_t end, 
                                             std::atomic<size_t>& processed_bytes, 
                                             size_t file_size, bool verbose, int thread_id) {
    ClusteringParseResult result;
    const char* ptr = data + begin;
    const char* chunk_end = data + end;
    
    // Adjust begin to start at beginning of a line if not at the start of the file
    if (begin > 0) {
        while (ptr < chunk_end && *(ptr-1) != '\n') {
            ptr++;
        }
    }
    
    // Process each line in the chunk
    while (ptr < chunk_end) {
        // Parse node ID
        uint64_t original_node_id = 0;
        while (ptr < chunk_end && *ptr >= '0' && *ptr <= '9') {
            original_node_id = original_node_id * 10 + (*ptr - '0');
            ptr++;
        }
        
        // Skip tab
        if (ptr < chunk_end && *ptr == '\t') ptr++;
        
        // Parse cluster ID as string
        std::string cluster_id;
        while (ptr < chunk_end && *ptr != '\n' && *ptr != '\r') {
            cluster_id.push_back(*ptr);
            ptr++;
        }
        
        // Skip to next line
        while (ptr < chunk_end && *ptr != '\n') ptr++;
        if (ptr < chunk_end) ptr++; // Skip newline
        
        // Add to result
        result.cluster_original_nodes[cluster_id].insert(original_node_id);
    }
    
    // Update progress
    if (verbose && thread_id == 0) {
        size_t current = processed_bytes.fetch_add(ptr - (data + begin));
        double progress = 100.0 * current / file_size;
        std::cout << "\rParsing clustering file: " << std::fixed << std::setprecision(1) 
                  << progress << "%" << std::flush;
    } else {
        processed_bytes.fetch_add(ptr - (data + begin));
    }
    
    return result;
}

// Load clustering from a TSV file (node_id, cluster_id)
Clustering load_clustering(const std::string& filename, const Graph& graph, 
                           int num_threads = std::thread::hardware_concurrency(),
                           bool verbose = false) {
    Clustering clustering;
    clustering.reset(graph.num_nodes);
    
    MappedFile file;
    if (!file.open(filename)) {
        std::cerr << "Failed to open clustering file: " << filename << std::endl;
        return clustering;
    }
    
    const char* data = file.data();
    size_t file_size = file.size();
    
    if (verbose) {
        std::cout << "Loading clustering from: " << filename << std::endl;
        std::cout << "Clustering file size: " << file_size / 1024 << " KB" << std::endl;
        std::cout << "Using " << num_threads << " threads for parsing" << std::endl;
    }
    
    // Split file into chunks for parallel processing
    std::vector<std::thread> threads;
    std::vector<ClusteringParseResult> chunk_results(num_threads);
    std::atomic<size_t> processed_bytes(0);
    
    size_t chunk_size = file_size / num_threads;
    if (chunk_size == 0) chunk_size = file_size;
    
    auto start_time = std::chrono::steady_clock::now();
    
    // Launch threads to parse the file
    for (int i = 0; i < num_threads; i++) {
        size_t begin = i * chunk_size;
        size_t end = (i == num_threads - 1) ? file_size : (i + 1) * chunk_size;
        
        threads.emplace_back([&, i, begin, end]() {
            chunk_results[i] = parse_clustering_chunk(data, begin, end, 
                                                       processed_bytes, file_size, 
                                                       verbose, i);
        });
    }
    
    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }
    
    if (verbose) {
        std::cout << std::endl; // End progress line
        
        auto parse_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_time);
        std::cout << "Parsing completed in " << parse_time.count() / 1000.0 << " seconds" << std::endl;
    }
    
    // Merge results from all threads
    std::unordered_map<std::string, std::unordered_set<uint64_t>> all_cluster_nodes;
    
    for (const auto& result : chunk_results) {
        for (const auto& [cluster_id, nodes] : result.cluster_original_nodes) {
            all_cluster_nodes[cluster_id].insert(nodes.begin(), nodes.end());
        }
    }
    
    if (verbose) {
        std::cout << "Found " << all_cluster_nodes.size() << " unique clusters" << std::endl;
    }
    
    // Convert original node IDs to internal node IDs and build the clustering
    size_t nodes_processed = 0;
    size_t nodes_found = 0;
    size_t total_nodes = 0;
    
    for (const auto& [cluster_id, original_nodes] : all_cluster_nodes) {
        total_nodes += original_nodes.size();
    }
    
    auto build_start_time = std::chrono::steady_clock::now();
    
    for (const auto& [cluster_id, original_nodes] : all_cluster_nodes) {
        // Add nodes to cluster
        for (uint64_t original_node_id : original_nodes) {
            auto it = graph.node_map.find(original_node_id);
            if (it != graph.node_map.end()) {
                uint32_t node_id = it->second;
                
                // Assign node to cluster using the new method
                clustering.assign_node_to_cluster(node_id, cluster_id);
                nodes_found++;
            } else {
                if (verbose) {
                    std::cerr << "Warning: Node " << original_node_id 
                              << " not found in graph" << std::endl;
                }
                clustering.assign_missing_node_to_cluster(original_node_id, cluster_id);
            }
            
            // Progress reporting
            nodes_processed++;
            if (verbose && nodes_processed % 1000000 == 0) {
                double progress = 100.0 * nodes_processed / total_nodes;
                std::cout << "\rBuilding clustering: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
    }
    
    if (verbose) {
        std::cout << std::endl; // End progress line
        
        auto build_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - build_start_time);
        std::cout << "Building completed in " << build_time.count() / 1000.0 << " seconds" << std::endl;
        
        std::cout << "Loaded " << clustering.get_non_empty_cluster_count() << " non-empty clusters" << std::endl;
        std::cout << "Assigned " << clustering.get_clustered_node_count() 
                  << " out of " << graph.num_nodes << " nodes to clusters" << std::endl;
        
        // Check for nodes not found in the graph
        std::cout << "Nodes found in graph: " << nodes_found 
                  << " out of " << total_nodes << " nodes in clustering file" << std::endl;
    }
    
    return clustering;
}

// Save a filtered clustering containing only nodes that exist in the given graph
bool save_filtered_clustering(const std::string& filename, 
    const Clustering& clustering,
    const Graph& graph,
    const std::unordered_map<uint32_t, uint32_t>& node_map,
    bool verbose = false) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return false;
    }

    if (verbose) {
        std::cout << "Saving filtered clustering to: " << filename << std::endl;
    }

    size_t nodes_written = 0;
    size_t nodes_skipped = 0;

    // For each node in the subgraph, write its cluster assignment if it exists
    for (uint32_t new_node_id = 0; new_node_id < graph.num_nodes; ++new_node_id) {
        // Get the original node ID from the graph's ID map
        uint64_t original_node_id = graph.id_map[new_node_id];

        // Find the original internal node ID from the node map
        uint32_t original_internal_id = UINT32_MAX;
        for (const auto& [orig_id, new_id] : node_map) {
            if (new_id == new_node_id) {
                original_internal_id = orig_id;
                break;
            }
        }

        if (original_internal_id == UINT32_MAX) {
            // Couldn't find the original internal ID - should not happen
            if (verbose) {
                std::cerr << "Warning: Could not find original internal ID for node " << new_node_id << std::endl;
            }
            nodes_skipped++;
            continue;
        }

        // Check if this node has a cluster assignment
        if (original_internal_id < clustering.node_to_cluster_idx.size() && 
            clustering.node_to_cluster_idx[original_internal_id] != UINT32_MAX) {

            uint32_t cluster_idx = clustering.node_to_cluster_idx[original_internal_id];
            const std::string& cluster_id = clustering.get_cluster_id(cluster_idx);

            // Write to file: original_node_id (from graph) and cluster_id
            outfile << original_node_id << "\t" << cluster_id << "\n";
            nodes_written++;
        } else {
            nodes_skipped++;
        }
    }

    outfile.close();

    if (verbose) {
        std::cout << "Written " << nodes_written << " node-cluster assignments" << std::endl;
        std::cout << "Skipped " << nodes_skipped << " nodes without cluster assignments" << std::endl;
    }

    return true;
}

#endif // CLUSTER_IO_H
