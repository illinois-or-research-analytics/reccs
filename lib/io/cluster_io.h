#ifndef CLUSTER_IO_H
#define CLUSTER_IO_H

#include <string>
#include <vector>
#include <atomic>
#include <mutex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include "../data_structures/clustering.h"
#include "../data_structures/graph.h"
#include "mapped_file.h"

// Helper struct for parallel loading (simplified - no mutex needed)
struct ClusteringParseResult {
    std::unordered_map<std::string, std::unordered_set<uint64_t>> cluster_original_nodes;
};

// Parse a range of the file - returns results without mutex issues
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
        // Skip empty lines and whitespace
        while (ptr < chunk_end && (*ptr == ' ' || *ptr == '\t' || *ptr == '\r' || *ptr == '\n')) {
            ptr++;
        }
        
        if (ptr >= chunk_end) break;
        
        // Parse node ID - check if we actually have digits
        const char* node_start = ptr;
        uint64_t original_node_id = 0;
        bool has_digits = false;
        
        while (ptr < chunk_end && *ptr >= '0' && *ptr <= '9') {
            original_node_id = original_node_id * 10 + (*ptr - '0');
            ptr++;
            has_digits = true;
        }
        
        // Validation: Must have parsed at least one digit and not overflow
        if (!has_digits) {
            // Skip this malformed line
            while (ptr < chunk_end && *ptr != '\n') ptr++;
            if (ptr < chunk_end) ptr++;
            continue;
        }
        
        // Check for potential overflow (node_id became smaller after adding digits)
        if (ptr - node_start > 19) { // uint64_t max is ~20 digits
            // Potential overflow, skip this line
            while (ptr < chunk_end && *ptr != '\n') ptr++;
            if (ptr < chunk_end) ptr++;
            continue;
        }
        
        // Must have a tab after node ID
        if (ptr >= chunk_end || *ptr != '\t') {
            // Malformed line, skip it
            while (ptr < chunk_end && *ptr != '\n') ptr++;
            if (ptr < chunk_end) ptr++;
            continue;
        }
        
        ptr++; // Skip tab
        
        // Parse cluster ID as string
        std::string cluster_id;
        while (ptr < chunk_end && *ptr != '\n' && *ptr != '\r' && *ptr != '\t') {
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
            while (ptr < chunk_end && *ptr != '\n') ptr++;
            if (ptr < chunk_end) ptr++;
            continue;
        }
        
        // Skip to end of line (handle both \r\n and \n)
        while (ptr < chunk_end && *ptr != '\n') ptr++;
        if (ptr < chunk_end) ptr++; // Skip newline
        
        // Validation: Check for reasonable node ID values
        if (original_node_id == 0) {
            if (verbose) {
                #pragma omp critical(warning_output)
                {
                    std::cerr << "Warning: Found node_id = 0 in cluster " << cluster_id << std::endl;
                }
            }
        }
        
        // Add to result (thread-local, no synchronization needed)
        result.cluster_original_nodes[cluster_id].insert(original_node_id);
    }
    
    // Update progress
    if (verbose && thread_id == 0) {
        size_t current = processed_bytes.fetch_add(ptr - (data + begin));
        double progress = 100.0 * current / file_size;
        #pragma omp critical(progress_output)
        {
            std::cout << "\rParsing clustering file: " << std::fixed << std::setprecision(1) 
                      << progress << "%" << std::flush;
        }
    } else {
        processed_bytes.fetch_add(ptr - (data + begin));
    }
    
    return result;
}

// Load clustering from a TSV file using OpenMP
Clustering load_clustering(const std::string& filename, const Graph& graph, 
                           int num_threads = 0, // 0 = use OpenMP default
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
    
    // Set number of threads
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    int actual_threads = omp_get_max_threads();
    
    if (verbose) {
        std::cout << "Loading clustering from: " << filename << std::endl;
        std::cout << "Clustering file size: " << file_size / 1024 << " KB" << std::endl;
        std::cout << "Using " << actual_threads << " threads for parsing" << std::endl;
    }
    
    // Shared data structures for parallel processing
    std::unordered_map<std::string, std::unordered_set<uint64_t>> all_cluster_nodes;
    std::mutex merge_mutex;
    std::atomic<size_t> processed_bytes(0);
    
    auto start_time = std::chrono::steady_clock::now();
    
    // Parse file in parallel using OpenMP
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int num_threads_actual = omp_get_num_threads();
        
        // Calculate chunk boundaries for this thread
        size_t chunk_size = file_size / num_threads_actual;
        size_t begin = thread_id * chunk_size;
        size_t end = (thread_id == num_threads_actual - 1) ? file_size : (thread_id + 1) * chunk_size;
        
        // Parse this thread's chunk
        ClusteringParseResult result = parse_clustering_chunk(data, begin, end, 
                                                               processed_bytes, file_size, 
                                                               verbose, thread_id);
        
        // Merge results into shared data structure
        {
            std::lock_guard<std::mutex> lock(merge_mutex);
            for (const auto& [cluster_id, nodes] : result.cluster_original_nodes) {
                all_cluster_nodes[cluster_id].insert(nodes.begin(), nodes.end());
            }
        }
    }
    
    if (verbose) {
        std::cout << std::endl; // End progress line
        
        auto parse_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start_time);
        std::cout << "Parsing completed in " << parse_time.count() / 1000.0 << " seconds" << std::endl;
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
    
    // Build clustering in parallel
    std::atomic<size_t> atomic_nodes_processed(0);
    std::atomic<size_t> atomic_nodes_found(0);
    std::mutex clustering_mutex;
    
    // Convert map to vector for easier parallel iteration
    std::vector<std::pair<std::string, std::unordered_set<uint64_t>>> cluster_vector;
    cluster_vector.reserve(all_cluster_nodes.size());
    for (const auto& pair : all_cluster_nodes) {
        cluster_vector.push_back(pair);
    }
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < cluster_vector.size(); ++i) {
        const auto& [cluster_id, original_nodes] = cluster_vector[i];
        
        size_t local_processed = 0;
        size_t local_found = 0;
        
        // Process nodes in this cluster
        for (uint64_t original_node_id : original_nodes) {
            auto it = graph.node_map.find(original_node_id);
            if (it != graph.node_map.end()) {
                uint32_t node_id = it->second;
                
                // Thread-safe assignment to clustering
                {
                    std::lock_guard<std::mutex> lock(clustering_mutex);
                    clustering.assign_node_to_cluster(node_id, cluster_id);
                }
                local_found++;
            } else {
                if (verbose) {
                    #pragma omp critical(node_warning)
                    {
                        std::cerr << "Warning: Node " << original_node_id 
                                  << " not found in graph" << std::endl;
                    }
                }
                {
                    std::lock_guard<std::mutex> lock(clustering_mutex);
                    clustering.assign_missing_node_to_cluster(original_node_id, cluster_id);
                }
            }
            
            local_processed++;
        }
        
        // Update atomic counters
        atomic_nodes_processed.fetch_add(local_processed);
        atomic_nodes_found.fetch_add(local_found);
        
        // Progress reporting (only from thread 0)
        if (verbose && omp_get_thread_num() == 0) {
            size_t current_processed = atomic_nodes_processed.load();
            if (current_processed % 1000000 == 0 || current_processed == total_nodes) {
                double progress = 100.0 * current_processed / total_nodes;
                #pragma omp critical(build_progress)
                {
                    std::cout << "\rBuilding clustering: " << std::fixed << std::setprecision(1) 
                              << progress << "%" << std::flush;
                }
            }
        }
    }
    
    nodes_processed = atomic_nodes_processed.load();
    nodes_found = atomic_nodes_found.load();
    
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

    std::atomic<size_t> nodes_written(0);
    std::atomic<size_t> nodes_skipped(0);
    std::mutex file_mutex;

    // Parallelize the writing process
    #pragma omp parallel for schedule(static, 1000)
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
                #pragma omp critical(save_warning)
                {
                    std::cerr << "Warning: Could not find original internal ID for node " << new_node_id << std::endl;
                }
            }
            nodes_skipped.fetch_add(1);
            continue;
        }

        // Check if this node has a cluster assignment
        if (original_internal_id < clustering.node_to_cluster_idx.size() && 
            clustering.node_to_cluster_idx[original_internal_id] != UINT32_MAX) {

            uint32_t cluster_idx = clustering.node_to_cluster_idx[original_internal_id];
            const std::string& cluster_id = clustering.get_cluster_id(cluster_idx);

            // Write to file with thread safety
            {
                std::lock_guard<std::mutex> lock(file_mutex);
                outfile << original_node_id << "\t" << cluster_id << "\n";
            }
            nodes_written.fetch_add(1);
        } else {
            nodes_skipped.fetch_add(1);
        }
    }

    outfile.close();

    if (verbose) {
        std::cout << "Written " << nodes_written.load() << " node-cluster assignments" << std::endl;
        std::cout << "Skipped " << nodes_skipped.load() << " nodes without cluster assignments" << std::endl;
    }

    return true;
}

#endif // CLUSTER_IO_H
