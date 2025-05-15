#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include <iomanip>
#include <filesystem>
#include <thread>
#include <unordered_set>
#include "../lib/data_structures/graph.h"
#include "../lib/io/graph_io.h"
#include "../lib/utils/orchestrator.h"

namespace fs = std::filesystem;

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

    // Check if clustering file is provided
    if (cluster_filename.empty()) {
        std::cerr << "Error: Clustering file (-c) is required" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    // Create a temp directory for intermediate files
    if (verbose) {
        std::cout << "Creating temporary directory for intermediate files..." << std::endl;
    }
    std::string temp_dir = "temp";
    if (!fs::exists(temp_dir)) {
        fs::create_directories(temp_dir);
    }
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads);
    
    // Create and run the orchestrator
    auto start_time = std::chrono::high_resolution_clock::now();
    Orchestrator orchestrator(graph_filename, cluster_filename, temp_dir, verbose);
    int result = orchestrator.run();
    
    if (result != 0) {
        std::cerr << "Orchestrator failed with error code " << result << std::endl;
        return result;
    }
    
    // Load the clustered SBM graph and clustering after orchestrator completes
    if (verbose) {
        std::cout << "Loading clustered SBM graph and clustering..." << std::endl;
    }
    
    // Define paths to the SBM-generated files
    std::string clustered_sbm_graph_path = temp_dir + "/clustered_sbm/syn_sbm.tsv";
    std::string clustered_clusters_path = temp_dir + "/non_singleton_clusters.tsv";
    
    // Check if files exist
    if (!fs::exists(clustered_sbm_graph_path)) {
        std::cerr << "Error: Clustered SBM graph file not found at: " << clustered_sbm_graph_path << std::endl;
        return 1;
    }
    
    if (!fs::exists(clustered_clusters_path)) {
        std::cerr << "Error: Clustered clustering file not found at: " << clustered_clusters_path << std::endl;
        return 1;
    }

    // Load the clustered SBM graph with node mapping
    Graph clustered_sbm_graph = GraphIO::read_tsv(clustered_sbm_graph_path, verbose);
    
    // Load the clustered clustering and map to internal IDs
    std::vector<int> clustered_clustering = GraphIO::read_clustering_mapped(
        clustered_clusters_path, clustered_sbm_graph, verbose);

    if (verbose) {
        std::cout << "Successfully loaded clustered SBM graph and clustering." << std::endl;
        clustered_sbm_graph.print_stats();
        
        // Count vertices with valid cluster assignments
        int assigned_vertices = 0;
        std::unordered_set<int> unique_clusters;
        
        for (size_t i = 0; i < clustered_clustering.size(); i++) {
            if (clustered_clustering[i] != -1) {
                assigned_vertices++;
                unique_clusters.insert(clustered_clustering[i]);
            }
        }
        
        std::cout << "Clustering statistics:" << std::endl;
        std::cout << "  Assigned vertices: " << assigned_vertices << " out of " 
                  << clustered_sbm_graph.num_vertices() << std::endl;
        std::cout << "  Unique clusters: " << unique_clusters.size() << std::endl;
        
        // Print timing information
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "Total execution time: " << duration << " seconds" << std::endl;
    }
    
    // Demonstrate working with original vs. internal IDs
    if (verbose) {
        std::cout << "\nNode mapping example:" << std::endl;
        
        // Show a sample of the node mapping
        const auto& node_map = clustered_sbm_graph.node_mapping();
        const auto& rev_map = clustered_sbm_graph.reverse_mapping();
        
        int sample_size = std::min(5, static_cast<int>(clustered_sbm_graph.num_vertices()));
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
            std::vector<int> neighbors = clustered_sbm_graph.get_neighbors_original(original_id);
            
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
    
    // You can now work with the loaded graph and clustering data
    // For example, you could:
    // 1. Analyze cluster properties
    // 2. Calculate graph metrics
    // 3. Visualize the graph structure
    // 4. Compare the SBM-generated graph with the original input
    
    return 0;
}
