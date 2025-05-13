#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <omp.h>

// CSR representation for undirected graph
struct CSRGraph {
    std::vector<uint32_t> row_ptr;  // Offsets for each node's edge list
    std::vector<uint32_t> col_idx;  // Target nodes
    
    // Node ID mapping (if needed)
    std::unordered_map<uint64_t, uint32_t> node_map;
    std::vector<uint64_t> id_map;
    
    // Graph info
    size_t num_nodes = 0;
    size_t num_edges = 0; // This counts each undirected edge once
};

// Memory-mapped file reader
class MappedFile {
private:
    void* mapped_data = nullptr;
    size_t file_size = 0;
    int fd = -1;

public:
    ~MappedFile() {
        if (mapped_data) munmap(mapped_data, file_size);
        if (fd >= 0) close(fd);
    }
    
    bool open(const std::string& filename) {
        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd == -1) return false;
        
        struct stat sb;
        if (fstat(fd, &sb) == -1) {
            close(fd);
            fd = -1;
            return false;
        }
        
        file_size = sb.st_size;
        mapped_data = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
        
        if (mapped_data == MAP_FAILED) {
            close(fd);
            fd = -1;
            mapped_data = nullptr;
            return false;
        }
        
        // Advise the kernel that we'll access this sequentially
        madvise(mapped_data, file_size, MADV_SEQUENTIAL);
        return true;
    }
    
    const char* data() const { return static_cast<const char*>(mapped_data); }
    size_t size() const { return file_size; }
};

CSRGraph load_undirected_tsv_edgelist_parallel(const std::string& filename, int num_threads = std::thread::hardware_concurrency()) {
    CSRGraph graph;
    MappedFile file;
    
    if (!file.open(filename)) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return graph;
    }
    
    const char* data = file.data();
    size_t file_size = file.size();
    
    // First pass: Count unique nodes and collect edges
    std::cout << "Step 1: Parsing file and collecting edges..." << std::endl;
    
    std::vector<std::thread> threads;
    std::vector<std::unordered_map<uint64_t, uint32_t>> local_node_maps(num_threads);
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> local_edges(num_threads);
    
    std::atomic<size_t> next_chunk_start(0);
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
    
    // Merge node maps and assign IDs
    std::cout << "Step 2: Building node ID mapping..." << std::endl;
    std::unordered_map<uint64_t, uint32_t>& node_map = graph.node_map;
    
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
    
    std::cout << "Found " << graph.num_nodes << " nodes and " << graph.num_edges << " undirected edges." << std::endl;
    
    // Count degrees in parallel
    std::cout << "Step 3: Counting node degrees..." << std::endl;
    
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
    
    // Build row pointers (prefix sum)
    std::cout << "Step 4: Building row pointers..." << std::endl;
    graph.row_ptr.resize(graph.num_nodes + 1);
    graph.row_ptr[0] = 0;
    
    for (size_t i = 0; i < graph.num_nodes; i++) {
        graph.row_ptr[i + 1] = graph.row_ptr[i] + degree[i].load(std::memory_order_relaxed);
    }
    
    // Prepare for parallel CSR building
    std::cout << "Step 5: Building CSR structure in parallel..." << std::endl;
    
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
    
    auto fill_csr = [&](int thread_id) {
        for (const auto& edge : local_edges[thread_id]) {
            uint32_t src_id = node_map[edge.first];
            uint32_t dst_id = node_map[edge.second];
            
            // Add edge src -> dst
            uint32_t pos1 = graph.row_ptr[src_id] + offsets[src_id].fetch_add(1, std::memory_order_relaxed);
            graph.col_idx[pos1] = dst_id;
            
            // Add edge dst -> src (undirected)
            uint32_t pos2 = graph.row_ptr[dst_id] + offsets[dst_id].fetch_add(1, std::memory_order_relaxed);
            graph.col_idx[pos2] = src_id;
        }
    };
    
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(fill_csr, i);
    }
    
    for (auto& t : threads) {
        t.join();
    }
    
    // Free memory from edge lists
    local_edges.clear();
    
    std::cout << "Loaded undirected graph with " << graph.num_nodes << " nodes and " 
              << graph.num_edges << " edges (" << total_directed_edges << " directed edges in CSR)" << std::endl;
    
    return graph;
}

// Sort adjacency lists in parallel
void sort_adjacency_lists_parallel(CSRGraph& graph, int num_threads) {
    std::cout << "Sorting adjacency lists..." << std::endl;
    
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        
        if (end > start) {
            std::sort(&graph.col_idx[start], &graph.col_idx[end]);
        }
    }
}

// Remove self-loops and duplicate edges in parallel
void clean_graph_parallel(CSRGraph& graph, int num_threads) {
    std::cout << "Removing self-loops and duplicate edges..." << std::endl;
    
    // First, sort adjacency lists to place duplicates adjacent to each other
    sort_adjacency_lists_parallel(graph, num_threads);
    
    // Count unique edges for each vertex
    std::vector<uint32_t> unique_counts(graph.num_nodes, 0);
    
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        
        if (start == end) continue; // No edges
        
        uint32_t last = UINT32_MAX; // Initialize to an impossible node ID
        uint32_t count = 0;
        
        for (uint32_t j = start; j < end; j++) {
            uint32_t neighbor = graph.col_idx[j];
            
            // Skip self-loops and duplicates
            if (neighbor != i && neighbor != last) {
                count++;
                last = neighbor;
            }
        }
        
        unique_counts[i] = count;
    }
    
    // Compute new row pointers based on unique counts
    std::vector<uint32_t> new_row_ptr(graph.num_nodes + 1);
    new_row_ptr[0] = 0;
    
    for (size_t i = 0; i < graph.num_nodes; i++) {
        new_row_ptr[i + 1] = new_row_ptr[i] + unique_counts[i];
    }
    
    // Allocate new column indices
    std::vector<uint32_t> new_col_idx(new_row_ptr.back());
    
    // Fill new column indices
    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < graph.num_nodes; i++) {
        uint32_t start = graph.row_ptr[i];
        uint32_t end = graph.row_ptr[i + 1];
        uint32_t new_start = new_row_ptr[i];
        
        if (start == end) continue; // No edges
        
        uint32_t last = UINT32_MAX; // Initialize to an impossible node ID
        uint32_t pos = new_start;
        
        for (uint32_t j = start; j < end; j++) {
            uint32_t neighbor = graph.col_idx[j];
            
            // Skip self-loops and duplicates
            if (neighbor != i && neighbor != last) {
                new_col_idx[pos++] = neighbor;
                last = neighbor;
            }
        }
    }
    
    // Update the graph
    graph.col_idx = std::move(new_col_idx);
    graph.row_ptr = std::move(new_row_ptr);
}

// Simple test function to validate the graph
void test_graph(const CSRGraph& graph) {
    std::cout << "Testing graph..." << std::endl;
    
    size_t total_edges = 0;
    for (size_t i = 0; i < graph.num_nodes; i++) {
        total_edges += graph.row_ptr[i + 1] - graph.row_ptr[i];
    }
    
    std::cout << "Total edges in CSR: " << total_edges << std::endl;
    
    // Check if all edges are valid
    bool valid = true;
    for (size_t i = 0; i < graph.num_nodes; i++) {
        for (size_t j = graph.row_ptr[i]; j < graph.row_ptr[i + 1]; j++) {
            uint32_t neighbor = graph.col_idx[j];
            if (neighbor >= graph.num_nodes) {
                std::cout << "Invalid edge: " << i << " -> " << neighbor << std::endl;
                valid = false;
                break;
            }
            
            // Check if the reverse edge exists
            bool found = false;
            for (size_t k = graph.row_ptr[neighbor]; k < graph.row_ptr[neighbor + 1]; k++) {
                if (graph.col_idx[k] == i) {
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                std::cout << "Missing reverse edge: " << neighbor << " -> " << i << std::endl;
                valid = false;
                break;
            }
        }
        
        if (!valid) break;
    }
    
    if (valid) {
        std::cout << "Graph is valid: all edges have corresponding reverse edges" << std::endl;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <edgelist.tsv> [num_threads]" << std::endl;
        return 1;
    }
    
    int num_threads = (argc > 2) ? std::stoi(argv[2]) : std::thread::hardware_concurrency();
    std::cout << "Using " << num_threads << " threads" << std::endl;
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads);
    
    auto total_start_time = std::chrono::steady_clock::now();
    
    // Load the graph with the fully parallel implementation
    auto load_start_time = std::chrono::steady_clock::now();
    CSRGraph graph = load_undirected_tsv_edgelist_parallel(argv[1], num_threads);
    auto load_end_time = std::chrono::steady_clock::now();
    
    // Optional: Clean the graph (remove self-loops and duplicates)
    auto clean_start_time = std::chrono::steady_clock::now();
    clean_graph_parallel(graph, num_threads);
    auto clean_end_time = std::chrono::steady_clock::now();
    
    // Optional: Validate the graph
    // test_graph(graph);
    
    auto total_end_time = std::chrono::steady_clock::now();
    
    // Print timing information
    auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(load_end_time - load_start_time);
    auto clean_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clean_end_time - clean_start_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time);
    
    std::cout << "Loading time: " << load_duration.count() / 1000.0 << " seconds" << std::endl;
    std::cout << "Cleaning time: " << clean_duration.count() / 1000.0 << " seconds" << std::endl;
    std::cout << "Total time: " << total_duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}
