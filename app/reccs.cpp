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
#include "../lib/io/g_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/io/requirements_io.h"
#include "../lib/io/degseq_io.h"
#include "../lib/utils/orchestrator.h"
#include "../lib/utils/edge_extractor.h"
#include "../lib/utils/statics.h"
#include "../lib/algorithm/enforce_degree_conn.h"
#include "../lib/algorithm/enforce_mincut.h"
#include "../lib/algorithm/deg_seq_matching.h"

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
    std::cerr << "  -c <clusters.tsv> Load clusters from TSV file (required)" << std::endl;
    std::cerr << "  -t <num_threads>  Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                Verbose mode: print detailed progress information" << std::endl;
    std::cerr << "  -o <output_file>  Output file (default: 'output.tsv')" << std::endl;
    std::cerr << "  -h, --help        Show this help message and exit" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Normal mode specific options:" << std::endl;
    std::cerr << "  <edgelist.tsv>                   Input graph edgelist file" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Checkpoint mode specific options:" << std::endl;
    std::cerr << "  --checkpoint                     Enable checkpoint mode (skip orchestrator)" << std::endl;
    std::cerr << "  --clustered-sbm <path>          Path to clustered SBM graph file" << std::endl;
    std::cerr << "  --unclustered-sbm <path>        Path to unclustered SBM graph file" << std::endl;
    std::cerr << "  --requirements <path>           Path to requirements CSV file" << std::endl;
    std::cerr << "  --degseq <path>                 Path to degree sequence JSON file" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Examples:" << std::endl;
    std::cerr << "  Normal mode:" << std::endl;
    std::cerr << "    " << program_name << " graph.tsv -c clusters.tsv -t 8 -v -o output.tsv" << std::endl;
    std::cerr << std::endl;
    std::cerr << "  Checkpoint mode:" << std::endl;
    std::cerr << "    " << program_name << " --checkpoint -c clusters.tsv \\" << std::endl;
    std::cerr << "      --clustered-sbm clustered.tsv --unclustered-sbm unclustered.tsv \\" << std::endl;
    std::cerr << "      --requirements requirements.csv --degseq degseq.json -v -o output.tsv" << std::endl;
}

struct CheckpointArgs {
    std::string clustered_sbm_path;
    std::string unclustered_sbm_path;
    std::string requirements_path;
    std::string degseq_path;
    
    bool are_all_provided() const {
        return !clustered_sbm_path.empty() && 
               !unclustered_sbm_path.empty() && 
               !requirements_path.empty() && 
               !degseq_path.empty();
    }
    
    std::vector<std::string> get_missing_args() const {
        std::vector<std::string> missing;
        if (clustered_sbm_path.empty()) missing.push_back("--clustered-sbm");
        if (unclustered_sbm_path.empty()) missing.push_back("--unclustered-sbm");
        if (requirements_path.empty()) missing.push_back("--requirements");
        if (degseq_path.empty()) missing.push_back("--degseq");
        return missing;
    }
    
    bool files_exist() const {
        return fs::exists(clustered_sbm_path) && 
               fs::exists(unclustered_sbm_path) && 
               fs::exists(requirements_path) && 
               fs::exists(degseq_path);
    }
    
    std::vector<std::string> get_missing_files() const {
        std::vector<std::string> missing;
        if (!fs::exists(clustered_sbm_path)) missing.push_back(clustered_sbm_path);
        if (!fs::exists(unclustered_sbm_path)) missing.push_back(unclustered_sbm_path);
        if (!fs::exists(requirements_path)) missing.push_back(requirements_path);
        if (!fs::exists(degseq_path)) missing.push_back(degseq_path);
        return missing;
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string graph_filename;
    std::string cluster_filename;
    std::string output_file = "output.tsv";
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    bool checkpoint_mode = false;
    CheckpointArgs checkpoint_args;
    
    // Check if first argument is --checkpoint
    if (std::string(argv[1]) == "--checkpoint") {
        checkpoint_mode = true;
        
        // Parse checkpoint mode arguments
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-v") {
                verbose = true;
            } else if (arg == "-t" && i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            } else if (arg == "-o" && i + 1 < argc) {
                output_file = argv[++i];
            } else if (arg == "-c" && i + 1 < argc) {
                cluster_filename = argv[++i];
            } else if (arg == "--clustered-sbm" && i + 1 < argc) {
                checkpoint_args.clustered_sbm_path = argv[++i];
            } else if (arg == "--unclustered-sbm" && i + 1 < argc) {
                checkpoint_args.unclustered_sbm_path = argv[++i];
            } else if (arg == "--requirements" && i + 1 < argc) {
                checkpoint_args.requirements_path = argv[++i];
            } else if (arg == "--degseq" && i + 1 < argc) {
                checkpoint_args.degseq_path = argv[++i];
            } else if (arg == "-h" || arg == "--help") {
                print_usage(argv[0]);
                return 0;
            } else {
                std::cerr << "Unknown checkpoint option: " << arg << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        
        // Validate checkpoint arguments (including degseq now)
        if (checkpoint_args.clustered_sbm_path.empty() || 
            checkpoint_args.unclustered_sbm_path.empty() ||
            checkpoint_args.requirements_path.empty() ||
            checkpoint_args.degseq_path.empty()) {
            std::cerr << "Error: In checkpoint mode, clustered-sbm, unclustered-sbm, requirements, and degseq arguments are required." << std::endl;
            std::vector<std::string> missing;
            if (checkpoint_args.clustered_sbm_path.empty()) missing.push_back("--clustered-sbm");
            if (checkpoint_args.unclustered_sbm_path.empty()) missing.push_back("--unclustered-sbm");
            if (checkpoint_args.requirements_path.empty()) missing.push_back("--requirements");
            if (checkpoint_args.degseq_path.empty()) missing.push_back("--degseq");
            
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
        
        // Check if required checkpoint files exist
        std::vector<std::string> files_to_check = {
            checkpoint_args.clustered_sbm_path,
            checkpoint_args.unclustered_sbm_path,
            checkpoint_args.requirements_path,
            checkpoint_args.degseq_path
        };
        
        std::vector<std::string> missing_files;
        for (const auto& file : files_to_check) {
            if (!fs::exists(file)) {
                missing_files.push_back(file);
            }
        }
        
        if (!missing_files.empty()) {
            std::cerr << "Error: Some checkpoint files do not exist:" << std::endl;
            for (const auto& file : missing_files) {
                std::cerr << "  " << file << std::endl;
            }
            return 1;
        }
        
        // Check if cluster file exists
        if (!fs::exists(cluster_filename)) {
            std::cerr << "Error: Cluster file does not exist: " << cluster_filename << std::endl;
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
            } else if (arg == "-t" && i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            } else if (arg == "-c" && i + 1 < argc) {
                cluster_filename = argv[++i];
            } else if (arg == "-o" && i + 1 < argc) {
                output_file = argv[++i];
            } else if (arg == "-h" || arg == "--help") {
                print_usage(argv[0]);
                return 0;
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
    
    // Set OpenMP threads
    omp_set_num_threads(num_threads);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::string clustered_sbm_graph_path;
    std::string unclustered_sbm_graph_path;
    std::string requirements_filename;
    std::string degseq_filename;
    
    if (checkpoint_mode) {
        if (verbose) {
            std::cout << "Running in checkpoint mode - skipping orchestrator..." << std::endl;
        }
        
        // Use checkpoint file paths directly
        clustered_sbm_graph_path = checkpoint_args.clustered_sbm_path;
        unclustered_sbm_graph_path = checkpoint_args.unclustered_sbm_path;
        requirements_filename = checkpoint_args.requirements_path;
        degseq_filename = checkpoint_args.degseq_path;
        
    } else {
        // Normal mode - run orchestrator
        if (verbose) {
            std::cout << "Running orchestrator..." << std::endl;
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
        
        // Create and run the orchestrator
        Orchestrator orchestrator(graph_filename, cluster_filename, temp_dir, verbose);
        int result = orchestrator.run();
        
        if (result != 0) {
            std::cerr << "Orchestrator failed with error code " << result << std::endl;
            return result;
        }
        
        // Set paths to orchestrator-generated files
        clustered_sbm_graph_path = temp_dir + "/clustered_sbm/syn_sbm.tsv";
        unclustered_sbm_graph_path = temp_dir + "/unclustered_sbm/syn_sbm.tsv";
        requirements_filename = temp_dir + "/clustered_stats.csv";
        degseq_filename = temp_dir + "/reference_degree_sequence.json";
    }
    
    // Load the clustered SBM graph and clustering
    if (verbose) {
        std::cout << "Loading clustered SBM graph and clustering..." << std::endl;
    }
    
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
    Clustering clustering = load_clustering(cluster_filename, clustered_sbm_graph, verbose);

    if (verbose) {
        std::cout << "Successfully loaded clustering with " << clustering.size() << " clusters." << std::endl;
    }

    // Load the cluster requirements
    if (verbose) {
        std::cout << "Loading cluster requirements from: " << requirements_filename << std::endl;
    }
    
    ConnectivityRequirementsLoader requirements_loader;
    if (!requirements_loader.load_from_csv(requirements_filename, verbose)) {
        std::cerr << "Error: Failed to load cluster requirements from " << requirements_filename << std::endl;
        return 1;
    }

    // Load the reference degree sequence
    std::shared_ptr<const std::vector<uint32_t>> reference_degree_sequence;
    
    if (verbose) {
        std::cout << "Loading reference degree sequence from: " << degseq_filename << std::endl;
    }
    
    json degseq_json = load_degseq_json(degseq_filename);
    if (degseq_json.empty()) {
        std::cerr << "Error: Failed to load reference degree sequence from " << degseq_filename << std::endl;
        return 1;
    }
    
    // Convert JSON array to vector<uint32_t> and wrap in shared_ptr
    auto sequence = std::make_shared<const std::vector<uint32_t>>(
        degseq_json.get<std::vector<uint32_t>>());
    reference_degree_sequence = sequence;
    
    if (verbose) {
        std::cout << "Successfully loaded reference degree sequence with " 
                  << reference_degree_sequence->size() << " degrees." << std::endl;
    }

    // Print loaded requirements statistics
    if (verbose) {
        requirements_loader.print_statistics();
    }

    // Load the graph task queue
    GraphTaskQueue task_queue;

    // Set up the task functions (removed degree sequence matching)
    task_queue.set_task_functions(
        // Connectivity enforcement (handles both degree and connectivity)
        [](Graph& g, uint32_t min_degree) {
            enforce_degree_and_connectivity(g, min_degree);
        },
        
        // WCC stitching
        [](Graph& g, uint32_t min_degree) {
            enforce_mincut(g, min_degree);
        }
    );

    // Initialize queue without degree sequence parameter
    task_queue.initialize_queue(clustered_sbm_graph, clustering, requirements_loader);
    
    if (verbose) {
        std::cout << "Initialized task queue with " << task_queue.queue_size() << " tasks." << std::endl;
    }

    // Process all tasks
    task_queue.process_all_tasks();

    // Retrieve completed subgraphs
    auto completed_subgraphs = task_queue.get_completed_subgraphs();
    if (verbose) {
        std::cout << "Processed " << completed_subgraphs.size() << " subgraphs." << std::endl;
    }

    // Output the added edges to a TSV file
    std::string added_edges_path;
    if (checkpoint_mode) {
        // In checkpoint mode, create a temp file for added edges
        added_edges_path = "temp_added_edges_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + ".tsv";
    } else {
        // Use temp directory from orchestrator
        std::string temp_dir = fs::path(clustered_sbm_graph_path).parent_path().parent_path();
        added_edges_path = temp_dir + "/added_edges.tsv";
    }
    
    auto newly_added_edges = EdgeExtractor::find_newly_added_edges(
        clustered_sbm_graph, completed_subgraphs);
    EdgeExtractor::write_edges_to_tsv_no_header(newly_added_edges, added_edges_path);

    if (verbose) {
        std::cout << "Wrote newly added edges to: " << added_edges_path << std::endl;
    }

    // Concatenate the clustered and unclustered SBM graphs, and the added edges
    std::ofstream output_stream(output_file);

    // Just dump all files into output
    for (const std::string& file : {clustered_sbm_graph_path, unclustered_sbm_graph_path, added_edges_path}) {
        std::ifstream in(file);
        output_stream << in.rdbuf();
    }

    if (verbose) {
        std::cout << "Concatenated clustered and unclustered SBM graphs, and added edges into: " 
                  << output_file << std::endl;
    }

    // Clean up temporary added_edges file if in checkpoint mode
    if (checkpoint_mode && fs::exists(added_edges_path)) {
        fs::remove(added_edges_path);
    }

    if (verbose) {
        // Print timing information
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "Total execution time: " << duration << " seconds" << std::endl;
    }
    
    return 0;
}
