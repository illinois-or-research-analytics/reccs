#include <iostream>
#include <chrono>
#include <omp.h>
#include "../lib/data_structures/graph.h"
#include "../lib/io/graph_io.h"

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
