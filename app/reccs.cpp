#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include "../lib/data_structures/graph.h"
#include "../lib/data_structures/clustering.h"
#include "../lib/io/graph_io.h"
#include "../lib/io/cluster_io.h"

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <edgelist.tsv> [options]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -t <num_threads>  Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                Verbose mode: print detailed progress information" << std::endl;
    std::cerr << "  -c <clusters.tsv> Load clusters from TSV file" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string graph_filename;
    std::string cluster_filename;
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
        } else if (arg == "-c" && i + 1 < argc) {
            cluster_filename = argv[++i];
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
        std::cout << "Loading graph from file: " << graph_filename << std::endl;
    }
    
    auto total_start_time = std::chrono::steady_clock::now();
    
    // Load the graph
    auto load_start_time = std::chrono::steady_clock::now();
    CSRGraph graph = load_undirected_tsv_edgelist_parallel(graph_filename, num_threads, verbose);
    auto load_end_time = std::chrono::steady_clock::now();
    
    // Clean the graph
    auto clean_start_time = std::chrono::steady_clock::now();
    clean_graph_parallel(graph, num_threads, verbose);
    auto clean_end_time = std::chrono::steady_clock::now();
    
    // Load clustering if specified
    if (!cluster_filename.empty()) {
        auto cluster_start_time = std::chrono::steady_clock::now();
        Clustering clustering = load_clustering(cluster_filename, graph, num_threads, verbose);
        auto cluster_end_time = std::chrono::steady_clock::now();
        
        auto cluster_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            cluster_end_time - cluster_start_time);
            
        if (verbose) {
            std::cout << "Cluster loading time: " << cluster_duration.count() / 1000.0 << " seconds" << std::endl;
            
            // Verify the clustering with verbose output
            std::cout << "Verifying clustering..." << std::endl;
            bool valid = clustering.verify(verbose);
            std::cout << "Clustering validation: " << (valid ? "VALID" : "INVALID") << std::endl;
            
            // Sample some clusters
            size_t sample_count = std::min(size_t(5), clustering.cluster_ids.size());
            std::cout << "Sample of loaded clusters:" << std::endl;
            for (size_t i = 0; i < sample_count; ++i) {
                const std::string& cluster_id = clustering.cluster_ids[i];
                const auto& nodes = clustering.get_cluster_nodes(cluster_id);
                std::cout << "  Cluster '" << cluster_id << "' (index " << i << "): " 
                          << nodes.size() << " nodes" << std::endl;
                
                // Show a few node IDs from this cluster
                if (!nodes.empty()) {
                    size_t node_count = 0;
                    std::cout << "    Nodes: ";
                    for (uint32_t node_id : nodes) {
                        std::cout << node_id << " ";
                        if (++node_count >= 5) break;
                    }
                    std::cout << (nodes.size() > 5 ? "..." : "") << std::endl;
                    
                    // Verify a sample node's assignment
                    uint32_t sample_node = *nodes.begin();
                    std::cout << "    Node " << sample_node << " assigned to cluster: " 
                              << clustering.get_node_cluster(sample_node) << std::endl;
                }
            }
        }
    }
    
    auto total_end_time = std::chrono::steady_clock::now();
    
    // Print timing information
    auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(load_end_time - load_start_time);
    auto clean_duration = std::chrono::duration_cast<std::chrono::milliseconds>(clean_end_time - clean_start_time);
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - total_start_time);
    
    if (verbose) {
        std::cout << "Graph loading time: " << load_duration.count() / 1000.0 << " seconds" << std::endl;
        std::cout << "Graph cleaning time: " << clean_duration.count() / 1000.0 << " seconds" << std::endl;
    }
    std::cout << "Total time: " << total_duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}
