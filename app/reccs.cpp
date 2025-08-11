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
#include "../lib/data_structures/graph_task_queue_with_degrees.h" // Modified include
#include "../lib/io/g_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/io/requirements_io.h"
#include "../lib/utils/orchestrator.h"
#include "../lib/utils/edge_extractor.h"
#include "../lib/utils/statics.h"

#include "../lib/algorithm/enforce_min_degree_with_budget.h"
#include "../lib/algorithm/enforce_connectivity_with_budget.h"
#include "../lib/algorithm/enforce_mincut_with_budget.h"

#include "../lib/algorithm/deg_seq_matching_v2.h"
#include "../lib/algorithm/deg_seq_matching_v1.h"

namespace fs = std::filesystem;
using json = nlohmann::json;

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <edgelist.tsv> [options]" << std::endl;
    std::cerr << "       " << program_name << " --checkpoint [checkpoint_options]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Note: In normal mode, <edgelist.tsv> must be the first argument." << std::endl;
    std::cerr << "      In checkpoint mode, --checkpoint must be the first argument." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Common options:" << std::endl;
    std::cerr << "  -c <clusters.tsv>   Load clusters from TSV file (required)" << std::endl;
    std::cerr << "  -t <num_threads>    Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                  Verbose mode: print detailed progress information" << std::endl;
    std::cerr << "  -o <output_file>    Output file (default: 'output.tsv')" << std::endl;
    std::cerr << "  -h, --help          Show this help message and exit" << std::endl;
    std::cerr << "  --v2                Use V2 degree sequence fitting with SBM." << std::endl;
    std::cerr << "  --cleanup           Clean up temporary files after execution" << std::endl;
    std::cerr << "  --tempname <name>   Use a custom temporary directory name (default: 'temp{timestamp}')" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Normal mode specific options:" << std::endl;
    std::cerr << "  <edgelist.tsv>                   Input graph edgelist file (used as empirical graph)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Checkpoint mode specific options:" << std::endl;
    std::cerr << "  --checkpoint                    Enable checkpoint mode (skip orchestrator)" << std::endl;
    std::cerr << "  --clustered-subgraph <path>     Path to the empirical clustered subgraph file" << std::endl;
    std::cerr << "  --clustered-sbm <path>          Path to clustered SBM graph file" << std::endl;
    std::cerr << "  --unclustered-sbm <path>        Path to unclustered SBM graph file" << std::endl;
    std::cerr << "  --requirements <path>           Path to requirements CSV file" << std::endl;
    std::cerr << "  --deficits <path>               Path to degree deficits JSON file (default: 'degree_deficits.json')" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "    " << program_name << " graph.tsv -c clusters.tsv -t 8 -v -o output.tsv" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Checkpoint mode with separate empirical graph:" << std::endl;
    std::cerr << "    " << program_name << " --checkpoint -c clusters.tsv -e empirical.tsv \\" << std::endl;
    std::cerr << "      --clustered-sbm clustered.tsv --unclustered-sbm unclustered.tsv \\" << std::endl;
    std::cerr << "      --requirements requirements.csv --deficits degree_deficits.json -v -o output.tsv" << std::endl;
    std::cerr << std::endl;
}

struct CheckpointArgs {
    std::string clustered_subgraph_path;
    std::string clustered_sbm_path;
    std::string unclustered_sbm_path;
    std::string requirements_path;
    std::string deficits_path;
    
    bool are_all_provided() const {
        return !clustered_subgraph_path.empty() &&
               !clustered_sbm_path.empty() && 
               !unclustered_sbm_path.empty() && 
               !requirements_path.empty() &&
               !deficits_path.empty();
    }
    
    std::vector<std::string> get_missing_args() const {
        std::vector<std::string> missing;
        if (clustered_subgraph_path.empty()) missing.push_back("--clustered-subgraph");
        if (clustered_sbm_path.empty()) missing.push_back("--clustered-sbm");
        if (unclustered_sbm_path.empty()) missing.push_back("--unclustered-sbm");
        if (requirements_path.empty()) missing.push_back("--requirements");
        if (deficits_path.empty()) missing.push_back("--deficits");
        return missing;
    }
    
    bool files_exist() const {
        return fs::exists(clustered_subgraph_path) &&
               fs::exists(clustered_sbm_path) && 
               fs::exists(unclustered_sbm_path) && 
               fs::exists(requirements_path);
    }

    std::vector<std::string> get_missing_files() const {
        std::vector<std::string> missing;
        if (!fs::exists(clustered_subgraph_path)) missing.push_back(clustered_subgraph_path);
        if (!fs::exists(clustered_sbm_path)) missing.push_back(clustered_sbm_path);
        if (!fs::exists(unclustered_sbm_path)) missing.push_back(unclustered_sbm_path);
        if (!fs::exists(requirements_path)) missing.push_back(requirements_path);
        if (!fs::exists(deficits_path)) missing.push_back(deficits_path);
        return missing;
    }
};

// Helper function to create node mappings
std::unordered_map<uint64_t, uint32_t> create_node_mapping(const Graph& graph) {
    std::unordered_map<uint64_t, uint32_t> mapping;
    for (uint32_t i = 0; i < graph.id_map.size(); ++i) {
        mapping[graph.id_map[i]] = i;
    }
    return mapping;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string graph_filename;
    std::string empirical_graph_filename; // Separate empirical graph
    std::string cluster_filename;
    std::string output_file = "output.tsv";
    std::string temp_dir;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    bool checkpoint_mode = false;
    bool use_v2 = false;
    bool cleanup = false; // Cleanup temporary files after execution
    bool tempname_provided = false;
    CheckpointArgs checkpoint_args;
    
    // Check if first argument is --checkpoint
    if (std::string(argv[1]) == "--checkpoint") {
        checkpoint_mode = true;
        
        // Parse checkpoint mode arguments
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-v") {
                verbose = true;
            } else if (arg == "--v2") {
                use_v2 = true;
            } else if (arg == "-t" && i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            } else if (arg == "-o" && i + 1 < argc) {
                output_file = argv[++i];
            } else if (arg == "-c" && i + 1 < argc) {
                cluster_filename = argv[++i];
            } else if (arg == "-e" && i + 1 < argc) {
                empirical_graph_filename = argv[++i];
            } else if (arg == "--clustered-subgraph" && i + 1 < argc) {
                checkpoint_args.clustered_subgraph_path = argv[++i];
            } else if (arg == "--clustered-sbm" && i + 1 < argc) {
                checkpoint_args.clustered_sbm_path = argv[++i];
            } else if (arg == "--unclustered-sbm" && i + 1 < argc) {
                checkpoint_args.unclustered_sbm_path = argv[++i];
            } else if (arg == "--requirements" && i + 1 < argc) {
                checkpoint_args.requirements_path = argv[++i];
            } else if (arg == "--deficits" && i + 1 < argc) {
                checkpoint_args.deficits_path = argv[++i];
            } else if (arg == "-h" || arg == "--help") {
                print_usage(argv[0]);
                return 0;
            } else if (arg == "--cleanup") {
                cleanup = true;
            } else if (arg == "--tempname" && i + 1 < argc) {
                tempname_provided = true;
                temp_dir = argv[++i];
            } else {
                std::cerr << "Unknown checkpoint option: " << arg << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        
        // Validate checkpoint arguments
        if (!checkpoint_args.are_all_provided()) {
            std::cerr << "Error: In checkpoint mode, all checkpoint arguments are required." << std::endl;
            auto missing = checkpoint_args.get_missing_args();
            std::cerr << "Missing arguments: ";
            for (size_t i = 0; i < missing.size(); ++i) {
                std::cerr << missing[i];
                if (i < missing.size() - 1) std::cerr << ", ";
            }
            std::cerr << std::endl;
            print_usage(argv[0]);
            return 1;
        }
        
        // Check if cluster file is provided
        if (cluster_filename.empty()) {
            std::cerr << "Error: Clustering file (-c) is required in checkpoint mode" << std::endl;
            print_usage(argv[0]);
            return 1;
        }
        
        // Check if files exist
        if (!checkpoint_args.files_exist()) {
            std::cerr << "Error: Some checkpoint files do not exist:" << std::endl;
            auto missing = checkpoint_args.get_missing_files();
            for (const auto& file : missing) {
                std::cerr << "  " << file << std::endl;
            }
            return 1;
        }
        
        if (!fs::exists(cluster_filename)) {
            std::cerr << "Error: Cluster file does not exist: " << cluster_filename << std::endl;
            return 1;
        }
        
        // If empirical graph specified, check it exists
        if (!empirical_graph_filename.empty() && !fs::exists(empirical_graph_filename)) {
            std::cerr << "Error: Empirical graph file does not exist: " << empirical_graph_filename << std::endl;
            return 1;
        }
        
    } else {
        // Normal mode - first argument is the graph filename
        graph_filename = argv[1];
        
        // Parse normal mode arguments
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-v") {
                verbose = true;
            } else if (arg == "--v2") {
                use_v2 = true;
            } else if (arg == "-t" && i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            } else if (arg == "-c" && i + 1 < argc) {
                cluster_filename = argv[++i];
            } else if (arg == "-e" && i + 1 < argc) {
                empirical_graph_filename = argv[++i];
            } else if (arg == "-o" && i + 1 < argc) {
                output_file = argv[++i];
            } else if (arg == "-h" || arg == "--help") {
                print_usage(argv[0]);
                return 0;
            } else if (arg == "--cleanup") {
                cleanup = true;
            } else if (arg == "--tempname" && i + 1 < argc) {
                tempname_provided = true;
                temp_dir = argv[++i];
            } else {
                std::cerr << "Unknown option: " << arg << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }

        // Check if clustering file is provided in normal mode
        if (cluster_filename.empty()) {
            std::cerr << "Error: Clustering file (-c) is required in normal mode" << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // In normal mode, if no separate empirical graph specified, use the input graph
    if (!checkpoint_mode && empirical_graph_filename.empty()) {
        empirical_graph_filename = graph_filename;
    }
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::string clustered_sbm_graph_path;
    std::string unclustered_sbm_graph_path;
    std::string requirements_filename;

    if (verbose) {
        std::cout << "Creating temporary directory for intermediate files..." << std::endl;
    }
    if (!tempname_provided) {
        temp_dir = "temp" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count());
    }

    if (checkpoint_mode) {
        if (verbose) {
            std::cout << "Running in checkpoint mode - skipping orchestrator..." << std::endl;
        }
        
        clustered_sbm_graph_path = checkpoint_args.clustered_sbm_path;
        unclustered_sbm_graph_path = checkpoint_args.unclustered_sbm_path;
        requirements_filename = checkpoint_args.requirements_path;
    } else {
        // Normal mode - run orchestrator
        if (verbose) {
            std::cout << "Running orchestrator..." << std::endl;
        }

        if (fs::exists(temp_dir)) {
            fs::remove_all(temp_dir);
        }

        fs::create_directories(temp_dir);
        
        Orchestrator orchestrator(graph_filename, cluster_filename, temp_dir, verbose);
        int result = orchestrator.run();
        
        if (result != 0) {
            std::cerr << "Orchestrator failed with error code " << result << std::endl;
            return result;
        }
        
        clustered_sbm_graph_path = temp_dir + "/clustered_sbm/syn_sbm.tsv";
        unclustered_sbm_graph_path = temp_dir + "/unclustered_sbm/syn_sbm.tsv";
        requirements_filename = temp_dir + "/clustered_stats.csv";
    }
    
    // Load the clustered SBM graph
    if (verbose) {
        std::cout << "Loading clustered SBM graph and clustering..." << std::endl;
    }
    
    if (!fs::exists(clustered_sbm_graph_path)) {
        std::cerr << "Error: Clustered SBM graph file not found at: " << clustered_sbm_graph_path << std::endl;
        return 1;
    }
    
    Graph clustered_sbm_graph = load_undirected_tsv_edgelist_parallel(
        clustered_sbm_graph_path, 1, verbose);

    omp_set_num_threads(num_threads);

    if (verbose) {
        std::cout << "Successfully loaded clustered SBM graph." << std::endl;
    }

    // Load the clustering
    Clustering clustering = load_clustering(cluster_filename, clustered_sbm_graph, verbose);

    if (verbose) {
        std::cout << "Successfully loaded clustering with " << clustering.size() << " clusters." << std::endl;
    }

    // Load requirements and degree sequence
    if (verbose) {
        std::cout << "Loading cluster requirements from: " << requirements_filename << std::endl;
    }
    
    ConnectivityRequirementsLoader requirements_loader;
    if (!requirements_loader.load_from_csv(requirements_filename, verbose)) {
        std::cerr << "Error: Failed to load cluster requirements from " << requirements_filename << std::endl;
        return 1;
    }

    if (verbose) {
        requirements_loader.print_statistics();
    }

    if (verbose) {
        std::cout << "\n=== STANDARD MODE TASK PROCESSING ===" << std::endl;
        std::cout << "Using degree-aware task queue with budget tracking..." << std::endl;
    }
    
    // Use degree-aware task queue
    GraphTaskQueueWithDegrees task_queue;

    // Initialize degree manager
    std::string deficits_path = checkpoint_mode ? checkpoint_args.deficits_path : temp_dir + "/degree_deficits.json";
    task_queue.initialize_degree_manager(deficits_path);

    // Set degree-aware task functions
    task_queue.set_task_functions(
        [](GraphTaskWithDegrees& task) {
            enforce_min_degree_with_budget(task);
        },
        [](GraphTaskWithDegrees& task) {
            enforce_connectivity_with_budget(task);
        },
        [](GraphTaskWithDegrees& task) {
            enforce_mincut_with_budget(task);
        }
    );

    task_queue.initialize_queue(clustered_sbm_graph, clustering, requirements_loader);
    
    if (verbose) {
        std::cout << "Initialized degree-aware task queue with " << task_queue.queue_size() << " tasks." << std::endl;
    }

    // Print initial degree budget statistics
    auto initial_stats = task_queue.get_degree_manager()->get_stats();
    if (verbose) {
        std::cout << "\n=== INITIAL DEGREE BUDGET STATISTICS ===" << std::endl;
        std::cout << "Available nodes: " << initial_stats.total_available_nodes << std::endl;
        std::cout << "Total degree budget: " << initial_stats.total_available_degrees << std::endl;
        std::cout << "Average budget per node: " << std::fixed << std::setprecision(2) 
                    << initial_stats.avg_available_degree << std::endl;
    }

    task_queue.process_all_tasks();

    // Print final degree budget statistics
    auto final_stats = task_queue.get_degree_manager()->get_stats();
    if (verbose) {
        std::cout << "\n=== FINAL DEGREE BUDGET STATISTICS ===" << std::endl;
        std::cout << "Remaining available nodes: " << final_stats.total_available_nodes << std::endl;
        std::cout << "Remaining degree budget: " << final_stats.total_available_degrees << std::endl;
        if (initial_stats.total_available_degrees > 0) {
            double utilization = 100.0 * (initial_stats.total_available_degrees - final_stats.total_available_degrees) / 
                                initial_stats.total_available_degrees;
            std::cout << "Budget utilization: " << std::fixed << std::setprecision(1) 
                        << utilization << "%" << std::endl;
        }
    }

    auto completed_subgraphs = task_queue.get_completed_subgraphs();
    if (verbose) {
        std::cout << "Processed " << completed_subgraphs.size() << " subgraphs in standard mode." << std::endl;
    }

    // Continue with edge extraction and final processing...
    auto newly_added_edges_raw = EdgeExtractor::find_newly_added_edges(
        clustered_sbm_graph, completed_subgraphs);

    std::vector<std::pair<uint32_t, uint32_t>> newly_added_edges = 
        EdgeExtractor::get_compressed_newly_added_edges(
            clustered_sbm_graph, newly_added_edges_raw, verbose);

    add_edges_batch(clustered_sbm_graph, newly_added_edges);

    // Modified clustered subgraph filename
    std::string temp_clustered_path = temp_dir + "/temp_clustered_modified.tsv";
    
    // Clustered subgraph path
    std::string clustered_subgraph_path = temp_dir + "/non_singleton_edges.tsv";
    std::string clustered_clusters_path = temp_dir + "/non_singleton_clusters.tsv";

    if (checkpoint_mode) {
        clustered_subgraph_path = checkpoint_args.clustered_subgraph_path;
        clustered_clusters_path = cluster_filename;
    }

    if (!use_v2) {
        if (verbose) {
            std::cout << "\nPerforming simple degree sequence matching (V1)..." << std::endl;
        }
        match_degree_sequence_v1(
            clustered_sbm_graph,
            clustered_subgraph_path,
            clustered_clusters_path,
            temp_dir + "/v1_output/",
            temp_clustered_path);
    } else {
        if (verbose) {
            std::cout << "\nPerforming SBM-based degree sequence matching (V2)..." << std::endl;
        }
        match_degree_sequence_v2(
            clustered_sbm_graph,
            clustered_subgraph_path,
            clustered_clusters_path, 
            temp_dir + "/sbm_output/",
            temp_clustered_path);
    }

    if (verbose) {
        std::cout << "Concatenating modified clustered and unclustered SBM graphs to: " << output_file << std::endl;
    }
    
    std::ofstream output_stream(output_file);
    
    std::ifstream clustered_in(temp_clustered_path);
    output_stream << clustered_in.rdbuf();
    clustered_in.close();
    
    std::ifstream unclustered_in(unclustered_sbm_graph_path);
    output_stream << unclustered_in.rdbuf();
    unclustered_in.close();
    
    output_stream.close();

    if (verbose) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "\n=== PROCESSING COMPLETE ===" << std::endl;
        std::cout << "Total execution time: " << duration << " seconds" << std::endl;
        std::cout << "Output written to: " << output_file << std::endl;
    }

    if (cleanup) {
        if (verbose) {
            std::cout << "Cleaning up temporary files..." << std::endl;
        }
        fs::remove_all(temp_dir);
    }
    
    return 0;
}
