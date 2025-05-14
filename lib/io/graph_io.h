#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>
#include <fstream>  
#include <iomanip>
#include <unordered_map>
#include <utility>
#include <chrono>
#include "../data_structures/graph.h"
#include "../data_structures/mapped_file.h"

DIGraph load_undirected_tsv_edgelist_parallel(const std::string& filename, int num_threads = std::thread::hardware_concurrency(), bool verbose = false) {
    DIGraph graph;
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
    
    // Initialize DI graph structure
    graph.init(graph.num_nodes);
    
    if (verbose) {
        std::cout << "Found " << graph.num_nodes << " nodes." << std::endl;
        std::cout << "Step 3: Building DI graph structure..." << std::endl;
    }
    
    // Count total edges
    size_t total_edges = 0;
    for (const auto& local_edge_list : local_edges) {
        total_edges += local_edge_list.size();
    }
    
    // Pre-allocate memory for edge arrays (2x for undirected)
    graph.src.reserve(total_edges * 2);
    graph.dst.reserve(total_edges * 2);
    
    // Fill edge arrays
    size_t edges_processed = 0;
    
    for (int thread_id = 0; thread_id < num_threads; thread_id++) {
        for (const auto& edge : local_edges[thread_id]) {
            uint32_t src_id = node_map[edge.first];
            uint32_t dst_id = node_map[edge.second];
            
            // Add edge in both directions (undirected graph)
            graph.src.push_back(src_id);
            graph.dst.push_back(dst_id);
            graph.src.push_back(dst_id);
            graph.dst.push_back(src_id);
            
            edges_processed++;
            
            if (verbose && edges_processed % 1000000 == 0) {
                double progress = static_cast<double>(edges_processed) / total_edges * 100.0;
                std::cout << "\rBuilding graph progress: " << std::fixed << std::setprecision(1) 
                          << progress << "%" << std::flush;
            }
        }
    }
    
    if (verbose) {
        std::cout << std::endl; // End the progress line
    }
    
    // Update edge count (counting undirected edges)
    graph.num_edges = total_edges;
    
    // Free memory from edge lists
    local_edges.clear();
    
    if (verbose) {
        std::cout << "Building DI graph index..." << std::endl;
    }
    
    // Build the vertex index
    graph.build_index(num_threads, verbose);
    
    if (verbose) {
        std::cout << "Loaded undirected graph with " << graph.num_nodes << " nodes and " 
                  << graph.num_edges << " edges (" << graph.src.size() << " directed edges in DI structure)" << std::endl;
    }
    
    return graph;
}

// Save graph to a TSV edgelist file
bool save_graph_edgelist(const std::string& filename, const DIGraph& graph, bool verbose = false) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return false;
    }
    
    if (verbose) {
        std::cout << "Saving graph to: " << filename << std::endl;
    }
    
    size_t edges_written = 0;
    size_t total_edges = graph.num_edges;
    
    // For each node, write its edges (but only in one direction to avoid duplicates)
    for (uint32_t src_id = 0; src_id < graph.num_nodes; ++src_id) {
        // Skip if node has no edges
        if (graph.str[src_id] == -1 || graph.nei[src_id] == 0) {
            continue;
        }
        
        uint64_t original_src_id = graph.id_map[src_id];
        
        // Get the starting position in the edge arrays
        uint32_t start_pos = graph.str[src_id];
        
        // Write out each edge (but only if src < dst to avoid duplicates)
        for (uint32_t i = 0; i < graph.nei[src_id]; ++i) {
            uint32_t dst_id = graph.dst[start_pos + i];
            
            // Only write edges where src_id < dst_id to avoid duplicates
            if (src_id < dst_id) {
                uint64_t original_dst_id = graph.id_map[dst_id];
                outfile << original_src_id << "\t" << original_dst_id << "\n";
                
                edges_written++;
                if (verbose && edges_written % 1000000 == 0) {
                    double progress = 100.0 * edges_written / total_edges;
                    std::cout << "\rSaving edges: " << std::fixed << std::setprecision(1) 
                              << progress << "%" << std::flush;
                }
            }
        }
    }
    
    if (verbose) {
        std::cout << std::endl; // End progress line
        std::cout << "Saved " << edges_written << " edges" << std::endl;
    }
    
    outfile.close();
    return true;
}


#endif // GRAPH_IO_H
