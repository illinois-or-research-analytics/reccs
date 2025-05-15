#include <iostream>
#include <chrono>
#include <string>
#include <thread>
#include <filesystem>
#include <unordered_set>
#include "../lib/data_structures/graph.h"
#include "../lib/io/graph_io.h"

namespace fs = std::filesystem;

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <edgelist.tsv> [options]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -t <num_threads>  Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                Verbose mode: print detailed information" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string graph_filename;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    
    // First argument is the graph filename
    graph_filename = argv[1];
    
    // Parse optional arguments
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-v") {
            verbose = true;
        } else if (arg == "-t" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    // Check if the graph file exists
    if (!fs::exists(graph_filename)) {
        std::cerr << "Error: Graph file not found: " << graph_filename << std::endl;
        return 1;
    }
    
    if (verbose) {
        std::cout << "Loading graph from: " << graph_filename << std::endl;
        std::cout << "Using " << num_threads << " threads" << std::endl;
    }
    
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Load the graph
    Graph graph = GraphIO::read_tsv(graph_filename, verbose);
    
    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    auto duration_sec = static_cast<double>(duration_ms) / 1000.0;
    
    // Print results
    std::cout << "Graph loaded successfully in " << duration_sec << " seconds." << std::endl;
    
    // Print graph statistics
    std::cout << "Graph statistics:" << std::endl;
    graph.print_stats();
    
    // If verbose, print additional information
    if (verbose) {
        std::cout << "\nNode mapping example:" << std::endl;
        
        // Show a sample of the node mapping
        const auto& node_map = graph.node_mapping();
        const auto& rev_map = graph.reverse_mapping();
        
        int sample_size = std::min(5, static_cast<int>(graph.num_vertices()));
        std::cout << "Sample of node mapping (original ID -> internal ID):" << std::endl;
        
        int counter = 0;
        for (const auto& pair : node_map) {
            if (counter++ >= sample_size) break;
            std::cout << "  " << pair.first << " -> " << pair.second << std::endl;
        }
        
        // Example of using original IDs to access the graph
        if (!node_map.empty()) {
            // Pick the first original ID from the mapping
            int original_id = rev_map[0];
            std::vector<int> neighbors = graph.get_neighbors_original(original_id);
            
            std::cout << "Neighbors of vertex " << original_id << " (original ID): ";
            if (neighbors.empty()) {
                std::cout << "none" << std::endl;
            } else {
                for (size_t i = 0; i < std::min(5UL, neighbors.size()); i++) {
                    std::cout << neighbors[i] << " ";
                }
                if (neighbors.size() > 5) {
                    std::cout << "... (" << neighbors.size() << " total)";
                }
                std::cout << std::endl;
            }
        }
    }
    
    return 0;
}
