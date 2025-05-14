#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <utility>
#include <chrono>
#include "../data_structures/graph.h"
#include "../data_structures/mapped_file.h"

CSRGraph load_undirected_tsv_edgelist_parallel(const std::string& filename, int num_threads = std::thread::hardware_concurrency(), bool verbose = false) {
    CSRGraph graph;
    MappedFile file;
    
    if (!file.open(filename)) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return graph;
    }
    
    const char* data = file.data();
    size_t file_size = file.size();
    
    if (verbose) {
        std::cout << "File size: " << file_size / (1024 * 1024) << " MB" << std::endl;
        std::cout << "Step 1: Parsing file and collecting edges..." << std::endl;
    }
    
    std::vector<std::thread> threads;
    std::vector<std::unordered_map<uint64_t, uint32_t>> local_node_maps(num_threads);
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> local_edges(num_threads);
    
    std::atomic<size_t> next_chunk_start(0);
    std::atomic<size_t> processed_bytes(0);
    size_t chunk_size = 64 * 1024 * 1024; // 64 MB chunks
    
    auto process_chunk = [&](int thread_id) {
        local_edges[thread_id].reserve(10000000); // Pre-allocate space for edges
        
        while (true) {
            // Get next chunk to process
            size_t chunk_begin = next_chunk_start.fetch_add(chunk_size);
            if (chunk_begin >= file_size) break;
            
            size_t chunk_end = std::min(chunk_begin + chunk_size, file_size);
            
            // Adjust chunk_begin to start at beginning of a line
            if (chunk_begin > 0) {
                while (chunk_begin < file_size && data[chunk_begin-1] != '\n') {
                    chunk_begin++;
                }
            }
            
            // Process the chunk
            const char* ptr = data + chunk_begin;
            const char* end = data + chunk_end;
            
            while (ptr < end) {
                // Parse source node
                uint64_t src = 0;
                while (ptr < end && *ptr >= '0' && *ptr <= '9') {
                    src = src * 10 + (*ptr - '0');
                    ptr++;
                }
                
                // Skip tab
                if (ptr < end && *ptr == '\t') ptr++;
                
                // Parse target node
                uint64_t dst = 0;
                while (ptr < end && *ptr >= '0' && *ptr <= '9') {
                    dst = dst * 10 + (*ptr - '0');
                    ptr++;
                }
                
                // Add to local maps
                local_node_maps[thread_id][src] = 0; // Temporary value
                local_node_maps[thread_id][dst] = 0; // Temporary value
                
                // Store the edge
                local_edges[thread_id].emplace_back(src, dst);
                
                // Skip to next line
                while (ptr < end && *ptr != '\n') ptr++;
                if (ptr < end) ptr++; // Skip newline
            }
            
            // Update progress
            processed_bytes.fetch_add(chunk_end - chunk_begin);
            
            if (verbose && thread_id == 0) {
                double progress = static_cast<double>(processed_bytes) / file_size * 100.0;
                std::cout << "\rParsing progress: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
    };
    
    // Launch threads for the first pass
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(process_chunk, i);
    }
    
    // Wait for threads to finish
    for (auto& t : threads) {
        t.join();
    }
    
    if (verbose) {
        std::cout << std::endl; // End the progress line
        std::cout << "Step 2: Building node ID mapping..." << std::endl;
    }
    
    std::unordered_map<uint64_t, uint32_t>& node_map = graph.node_map;
    
    // Report on thread distribution if verbose
    if (verbose) {
        for (int i = 0; i < num_threads; i++) {
            std::cout << "Thread " << i << " processed " 
                      << local_edges[i].size() << " edges and found "
                      << local_node_maps[i].size() << " unique nodes." << std::endl;
        }
    }
    
    for (const auto& local_map : local_node_maps) {
        for (const auto& entry : local_map) {
            node_map[entry.first] = 0;
        }
    }
    
    // Clear local node maps to free memory
    local_node_maps.clear();
    
    // Assign sequential IDs
    uint32_t next_id = 0;
    for (auto& entry : node_map) {
        entry.second = next_id++;
    }
    
    graph.num_nodes = node_map.size();
    
    // Create reverse mapping
    graph.id_map.resize(graph.num_nodes);
    for (const auto& entry : node_map) {
        graph.id_map[entry.second] = entry.first;
    }
    
    // Calculate total edges (undirected)
    for (const auto& local_edge_list : local_edges) {
        graph.num_edges += local_edge_list.size();
    }
    
    if (verbose) {
        std::cout << "Found " << graph.num_nodes << " nodes and " 
                  << graph.num_edges << " undirected edges." << std::endl;
        std::cout << "Step 3: Counting node degrees..." << std::endl;
    }
    
    // Use atomic vector for thread-safe degree counting
    std::vector<std::atomic<uint32_t>> degree(graph.num_nodes);
    for (auto& d : degree) {
        d.store(0, std::memory_order_relaxed);
    }
    
    // Launch threads to count degrees
    threads.clear();
    
    auto count_degrees = [&](int thread_id) {
        for (const auto& edge : local_edges[thread_id]) {
            uint32_t src_id = node_map[edge.first];
            uint32_t dst_id = node_map[edge.second];
            
            // Increment degree for both source and destination (undirected graph)
            degree[src_id].fetch_add(1, std::memory_order_relaxed);
            degree[dst_id].fetch_add(1, std::memory_order_relaxed);
        }
    };
    
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(count_degrees, i);
    }
    
    for (auto& t : threads) {
        t.join();
    }
    
    if (verbose) {
        std::cout << "Step 4: Building row pointers..." << std::endl;
    }
    
    // Build row pointers (prefix sum)
    graph.row_ptr.resize(graph.num_nodes + 1);
    graph.row_ptr[0] = 0;
    
    for (size_t i = 0; i < graph.num_nodes; i++) {
        graph.row_ptr[i + 1] = graph.row_ptr[i] + degree[i].load(std::memory_order_relaxed);
    }
    
    if (verbose) {
        std::cout << "Step 5: Building CSR structure in parallel..." << std::endl;
    }
    
    // Prepare for parallel CSR building
    
    // Allocate column indices vector
    size_t total_directed_edges = graph.row_ptr.back();
    graph.col_idx.resize(total_directed_edges);
    
    // Use atomic offsets for thread-safe insertion
    std::vector<std::atomic<uint32_t>> offsets(graph.num_nodes);
    for (size_t i = 0; i < graph.num_nodes; i++) {
        offsets[i].store(0, std::memory_order_relaxed);
    }
    
    // Launch threads to fill CSR
    threads.clear();
    std::atomic<size_t> edges_processed(0);
    
    auto fill_csr = [&](int thread_id) {
        size_t local_edges_count = local_edges[thread_id].size();
        size_t local_processed = 0;
        
        for (const auto& edge : local_edges[thread_id]) {
            uint32_t src_id = node_map[edge.first];
            uint32_t dst_id = node_map[edge.second];
            
            // Add edge src -> dst
            uint32_t pos1 = graph.row_ptr[src_id] + offsets[src_id].fetch_add(1, std::memory_order_relaxed);
            graph.col_idx[pos1] = dst_id;
            
            // Add edge dst -> src (undirected)
            uint32_t pos2 = graph.row_ptr[dst_id] + offsets[dst_id].fetch_add(1, std::memory_order_relaxed);
            graph.col_idx[pos2] = src_id;
            
            local_processed++;
            
            // Update progress
            if (verbose && thread_id == 0 && local_processed % 1000000 == 0) {
                size_t total_processed = edges_processed.fetch_add(1000000);
                double progress = static_cast<double>(total_processed) / graph.num_edges * 100.0;
                std::cout << "\rCSR building progress: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
        
        // Account for any remaining edges
        if (local_processed % 1000000 != 0) {
            edges_processed.fetch_add(local_processed % 1000000);
        }
    };
    
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(fill_csr, i);
    }
    
    for (auto& t : threads) {
        t.join();
    }
    
    if (verbose) {
        std::cout << std::endl; // End the progress line
    }
    
    // Free memory from edge lists
    local_edges.clear();
    
    if (verbose) {
        std::cout << "Loaded undirected graph with " << graph.num_nodes << " nodes and " 
                  << graph.num_edges << " edges (" << total_directed_edges << " directed edges in CSR)" << std::endl;
    }
    
    return graph;
}

#endif // GRAPH_IO_H
