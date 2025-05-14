#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include "../lib/data_structures/graph.h"
#include "../lib/io/graph_io.h"

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <edgelist.tsv> [options]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -t <num_threads>  Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                Verbose mode: print detailed progress information" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string filename;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    
    // First argument is the filename
    filename = argv[1];
    
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
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads);
    
    if (verbose) {
        std::cout << "Using " << num_threads << " threads" << std::endl;
        std::cout << "Loading graph from file: " << filename << std::endl;
    }
    
    auto total_start_time = std::chrono::steady_clock::now();
    
    // Load the graph with the fully parallel implementation
    auto load_start_time = std::chrono::steady_clock::now();
    CSRGraph graph = load_undirected_tsv_edgelist_parallel(filename, num_threads, verbose);
    auto load_end_time = std::chrono::steady_clock::now();
    
    // Optional: Clean the graph (remove self-loops and duplicates)
    auto clean_start_time = std::chrono::steady_clock::now();
    clean_graph_parallel(graph, num_threads, verbose);
    auto clean_end_time = std::chrono::steady_clock::now();

    // Get the total time taken
    auto total_end_time = std::chrono::steady_clock::now();
    
    // Print timing information
    auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(load_end_time - load_start_time);
    auto clean_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clean_end_time - clean_start_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time);
    
    if (verbose) {
        std::cout << "Loading time: " << load_duration.count() / 1000.0 << " seconds" << std::endl;
        std::cout << "Cleaning time: " << clean_duration.count() / 1000.0 << " seconds" << std::endl;
    }
    std::cout << "Total time: " << total_duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}
