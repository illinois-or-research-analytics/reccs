#ifndef CLUSTER_IO_H
#define CLUSTER_IO_H

#include <string>
#include <vector>
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

// Load clustering from a TSV file (single-threaded, efficient)
Clustering load_clustering(const std::string& filename, const Graph& graph, 
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
    }
    
    auto start_time = std::chrono::steady_clock::now();
    
    // Step 1: Parse file and collect cluster assignments
    std::unordered_map<std::string, std::unordered_set<uint64_t>> all_cluster_nodes;
    
    const char* ptr = data;
    const char* end = data + file_size;
    size_t lines_processed = 0;
    
    if (verbose) {
        std::cout << "Step 1: Parsing clustering file..." << std::endl;
    }
    
    while (ptr < end) {
        // Skip empty lines and whitespace
        while (ptr < end && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r' || *ptr == '\n')) {
            ptr++;
        }
        
        if (ptr >= end) break;
        
        // Parse node ID - check if we actually have digits
        const char* node_start = ptr;
        uint64_t original_node_id = 0;
        bool has_digits = false;
        
        while (ptr < end && *ptr >= '0' && *ptr <= '9') {
            original_node_id = original_node_id * 10 + (*ptr - '0');
            ptr++;
            has_digits = true;
        }
        
        // Validation: Must have parsed at least one digit and not overflow
        if (!has_digits) {
            // Skip this malformed line
            while (ptr < end && *ptr != '\n') ptr++;
            if (ptr < end) ptr++;
            continue;
        }
        
        // Check for potential overflow (node_id became smaller after adding digits)
        if (ptr - node_start > 19) { // uint64_t max is ~20 digits
            // Potential overflow, skip this line
            while (ptr < end && *ptr != '\n') ptr++;
            if (ptr < end) ptr++;
            continue;
        }
        
        // Must have a tab after node ID
        if (ptr >= end || *ptr != '\t') {
            // Malformed line, skip it
            while (ptr < end && *ptr != '\n') ptr++;
            if (ptr < end) ptr++;
            continue;
        }
        
        ptr++; // Skip tab
        
        // Parse cluster ID as string
        std::string cluster_id;
        cluster_id.reserve(32); // Reserve space for typical cluster ID
        while (ptr < end && *ptr != '\n' && *ptr != '\r' && *ptr != '\t') {
            cluster_id.push_back(*ptr);
            ptr++;
        }
        
        // Trim trailing whitespace from cluster_id
        while (!cluster_id.empty() && (cluster_id.back() == ' ' || cluster_id.back() == '\t')) {
            cluster_id.pop_back();
        }
        
        // Validation: cluster_id must not be empty
        if (cluster_id.empty()) {
            // Skip this malformed line
            while (ptr < end && *ptr != '\n') ptr++;
            if (ptr < end) ptr++;
            continue;
        }
        
        // Skip to end of line (handle both \r\n and \n)
        while (ptr < end && *ptr != '\n') ptr++;
        if (ptr < end) ptr++; // Skip newline
        
        // Validation: Check for reasonable node ID values
        if (original_node_id == 0 && verbose) {
            std::cerr << "Warning: Found node_id = 0 in cluster " << cluster_id << std::endl;
        }
        
        // Add to cluster assignments
        all_cluster_nodes[cluster_id].insert(original_node_id);
        
        lines_processed++;
        
        // Progress reporting
        if (verbose && lines_processed % 1000000 == 0) {
            size_t current_pos = ptr - data;
            double progress = 100.0 * current_pos / file_size;
            std::cout << "\rParsing progress: " << std::fixed << std::setprecision(1) 
                      << progress << "%" << std::flush;
        }
    }
    
    if (verbose) {
        std::cout << std::endl; // End progress line
        
        auto parse_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_time);
        std::cout << "Parsing completed in " << parse_time.count() / 1000.0 << " seconds" << std::endl;
        std::cout << "Processed " << lines_processed << " lines" << std::endl;
        std::cout << "Found " << all_cluster_nodes.size() << " unique clusters" << std::endl;
    }
    
    // Step 2: Convert original node IDs to internal node IDs and build the clustering
    size_t nodes_processed = 0;
    size_t nodes_found = 0;
    size_t total_nodes = 0;
    
    for (const auto& [cluster_id, original_nodes] : all_cluster_nodes) {
        total_nodes += original_nodes.size();
    }
    
    auto build_start_time = std::chrono::steady_clock::now();
    
    if (verbose) {
        std::cout << "Step 2: Building clustering structure..." << std::endl;
    }
    
    for (const auto& [cluster_id, original_nodes] : all_cluster_nodes) {
        for (uint64_t original_node_id : original_nodes) {
            auto it = graph.node_map.find(original_node_id);
            if (it != graph.node_map.end()) {
                uint32_t node_id = it->second;
                clustering.assign_node_to_cluster(node_id, cluster_id);
                nodes_found++;
            } else {
                if (verbose) {
                    std::cerr << "Warning: Node " << original_node_id 
                              << " not found in graph" << std::endl;
                }
                clustering.assign_missing_node_to_cluster(original_node_id, cluster_id);
            }
            
            nodes_processed++;
            
            // Progress reporting
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

        // Debug: if not all nodes were assigned, report missing nodes
        if (clustering.get_clustered_node_count() < graph.num_nodes) {
            std::cout << "Warning: Not all nodes were assigned to clusters!" << std::endl;
            std::cout << "List of missing nodes:" << std::endl;
            for (uint32_t node_id = 0; node_id < graph.num_nodes; ++node_id) {
                if (clustering.node_to_cluster_idx[node_id] == UINT32_MAX) {
                    std::cout << "Node " << graph.id_map[node_id] << "(" << node_id << ")"
                              << " not assigned to any cluster" << std::endl;
                }
            }
        }
        
        // Check for nodes not found in the graph
        std::cout << "Nodes found in graph: " << nodes_found 
                  << " out of " << total_nodes << " nodes in clustering file" << std::endl;
        
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_time);
        std::cout << "Total loading time: " << total_time.count() / 1000.0 << " seconds" << std::endl;
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

    // Create reverse lookup map for efficiency (only if needed)
    std::unordered_map<uint32_t, uint32_t> reverse_node_map;
    for (const auto& [orig_id, new_id] : node_map) {
        reverse_node_map[new_id] = orig_id;
    }

    for (uint32_t new_node_id = 0; new_node_id < graph.num_nodes; ++new_node_id) {
        // Get the original node ID from the graph's ID map
        uint64_t original_node_id = graph.id_map[new_node_id];

        // Find the original internal node ID from the reverse map
        auto reverse_it = reverse_node_map.find(new_node_id);
        if (reverse_it == reverse_node_map.end()) {
            // Couldn't find the original internal ID - should not happen
            if (verbose) {
                std::cerr << "Warning: Could not find original internal ID for node " << new_node_id << std::endl;
            }
            nodes_skipped++;
            continue;
        }
        
        uint32_t original_internal_id = reverse_it->second;

        // Check if this node has a cluster assignment
        if (original_internal_id < clustering.node_to_cluster_idx.size() && 
            clustering.node_to_cluster_idx[original_internal_id] != UINT32_MAX) {

            uint32_t cluster_idx = clustering.node_to_cluster_idx[original_internal_id];
            const std::string& cluster_id = clustering.get_cluster_id(cluster_idx);

            outfile << original_node_id << "\t" << cluster_id << "\n";
            nodes_written++;
        } else {
            nodes_skipped++;
        }
        
        // Progress reporting
        if (verbose && (new_node_id + 1) % 1000000 == 0) {
            double progress = 100.0 * (new_node_id + 1) / graph.num_nodes;
            std::cout << "\rSaving progress: " << std::fixed << std::setprecision(1) 
                      << progress << "%" << std::flush;
        }
    }

    outfile.close();

    if (verbose) {
        std::cout << std::endl; // End progress line
        std::cout << "Written " << nodes_written << " node-cluster assignments" << std::endl;
        std::cout << "Skipped " << nodes_skipped << " nodes without cluster assignments" << std::endl;
    }

    return true;
}

#endif // CLUSTER_IO_H
