#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include <iomanip>
#include <filesystem>
#include <thread>
#include <unordered_set>
#include "../lib/data_structures/graph.h"
#include "../lib/data_structures/clustering.h"
#include "../lib/data_structures/graph_task_queue.h"
#include "../lib/io/graph_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/io/requirements_io.h"
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
    std::string temp_dir = "temp" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count());

    // Remove existing temp directory if it exists
    if (fs::exists(temp_dir)) {
        fs::remove_all(temp_dir);
    }

    fs::create_directories(temp_dir);
    
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
    
    // Check if files exist
    if (!fs::exists(clustered_sbm_graph_path)) {
        std::cerr << "Error: Clustered SBM graph file not found at: " << clustered_sbm_graph_path << std::endl;
        return 1;
    }
    
    // Load the clustered SBM graph with node mapping
    Graph clustered_sbm_graph = load_undirected_tsv_edgelist_parallel(
        clustered_sbm_graph_path, num_threads, verbose);

    if (verbose) {
        std::cout << "Successfully loaded clustered SBM graph." << std::endl;
    }

    // Load the clustering from the specified file
    Clustering clustering = load_clustering(cluster_filename, clustered_sbm_graph, num_threads, verbose);

    if (verbose) {
        std::cout << "Successfully loaded clustering with " << clustering.size() << " clusters." << std::endl;
    }

    // Load the cluster requirements
    std::string requirements_filename = temp_dir + "/clustered_stats.csv";
    if (verbose) {
        std::cout << "Loading cluster requirements from: " << requirements_filename << std::endl;
    }
    
    ConnectivityRequirementsLoader requirements_loader;
    if (!requirements_loader.load_from_csv(requirements_filename, verbose)) {
        std::cerr << "Error: Failed to load cluster requirements from " << requirements_filename << std::endl;
        return 1;
    }

    // Print loaded requirements statistics
    if (verbose) {
        requirements_loader.print_statistics();
    }

    // Load the graph task queue
    GraphTaskQueue task_queue;
    task_queue.initialize_queue(clustered_sbm_graph, clustering);

    if (verbose) {
        std::cout << "Initialized task queue with " << task_queue.queue_size() << " tasks." << std::endl;
    }

    if (verbose) {
        // Print timing information
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "Total execution time: " << duration << " seconds" << std::endl;
    }
    
    return 0;
}
