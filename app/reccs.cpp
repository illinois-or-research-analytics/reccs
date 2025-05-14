#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include <iomanip>
#include <sys/wait.h>
#include <unistd.h>
#include <filesystem>
#include "../lib/data_structures/graph.h"
#include "../lib/data_structures/clustering.h"
#include "../lib/io/graph_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/algorithms/graph_splitter.h"
#include "../lib/algorithms/graph_merger.h"
#include "../lib/algorithms/sbm.h"
#include "../lib/algorithms/degree_enforcer.h"

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
    
    // Create a base output directory name from the input filename
    std::string base_name;
    size_t dot_pos = graph_filename.find_last_of('.');
    size_t slash_pos = graph_filename.find_last_of('/');
    
    if (dot_pos != std::string::npos && (slash_pos == std::string::npos || slash_pos < dot_pos)) {
        base_name = graph_filename.substr(0, dot_pos);
    } else {
        base_name = graph_filename;
    }
    
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
    
    // Create temporary directory for outputs
    std::string temp_dir = base_name + "_temp";
    if (!fs::exists(temp_dir)) {
        fs::create_directories(temp_dir);
    }
    
    if (verbose) {
        std::cout << "Created temporary directory: " << temp_dir << std::endl;
    }
    
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
                }
            }
        }
        
        // Split the graph
        auto split_start_time = std::chrono::steady_clock::now();
        
        auto [clustered_graph, unclustered_graph, clustered_node_map, unclustered_node_map] = 
            GraphSplitter::split_by_clustering(graph, clustering, verbose);
        
        auto split_end_time = std::chrono::steady_clock::now();
        auto split_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            split_end_time - split_start_time);
            
        if (verbose) {
            std::cout << "Graph splitting time: " << split_duration.count() / 1000.0 << " seconds" << std::endl;
            
            // Print stats about the split graphs
            std::cout << "Original graph: " << graph.num_nodes << " nodes, " 
                      << graph.num_edges << " edges" << std::endl;
            std::cout << "Clustered graph: " << clustered_graph.num_nodes << " nodes, " 
                      << clustered_graph.num_edges << " edges" << std::endl;
            std::cout << "Unclustered graph: " << unclustered_graph.num_nodes << " nodes, " 
                      << unclustered_graph.num_edges << " edges" << std::endl;
            
            // Compute some basic metrics
            double clustered_ratio = (double)clustered_graph.num_nodes / graph.num_nodes * 100.0;
            double clustered_edge_ratio = (double)clustered_graph.num_edges / graph.num_edges * 100.0;
            
            std::cout << "Percentage of nodes in non-singleton clusters: " 
                      << std::fixed << std::setprecision(2) << clustered_ratio << "%" << std::endl;
            std::cout << "Percentage of edges in non-singleton clusters: " 
                      << std::fixed << std::setprecision(2) << clustered_edge_ratio << "%" << std::endl;
        }
        
        // Save the split graphs to temporary directory
        std::string clustered_filename = temp_dir + "/clustered.tsv";
        save_graph_edgelist(clustered_filename, clustered_graph, verbose);
        
        std::string unclustered_filename = temp_dir + "/unclustered.tsv";
        save_graph_edgelist(unclustered_filename, unclustered_graph, verbose);
        
        if (verbose) {
            std::cout << "Saved clustered graph to: " << clustered_filename << std::endl;
            std::cout << "Saved unclustered graph to: " << unclustered_filename << std::endl;
        }
        
        // Generate SBM models for both subgraphs
        auto sbm_start_time = std::chrono::steady_clock::now();
        
        if (verbose) {
            std::cout << "Generating SBM models for clustered and unclustered subgraphs..." << std::endl;
        }
        
        // Create output directories for SBM models
        std::string clustered_sbm_dir = temp_dir + "/clustered_sbm";
        std::string unclustered_sbm_dir = temp_dir + "/unclustered_sbm";
        
        // Create filtered clustering files for each subgraph
        std::string clustered_clustering_file = temp_dir + "/clustered_clustering.tsv";
        std::string unclustered_clustering_file = temp_dir + "/unclustered_clustering.tsv";
        
        // Save filtered clusterings
        if (verbose) {
            std::cout << "Creating filtered clustering files for subgraphs..." << std::endl;
        }
        
        // Save filtered clustering files
        save_filtered_clustering(clustered_clustering_file, clustering, clustered_graph, clustered_node_map, verbose);
        save_filtered_clustering(unclustered_clustering_file, clustering, unclustered_graph, unclustered_node_map, verbose);
        
        // Fork processes to generate SBM models in parallel
        pid_t clustered_pid = SBMGenerator::fork_generate_sbm(
            clustered_graph, clustered_clustering_file, clustered_sbm_dir, "clustered", verbose);
            
        pid_t unclustered_pid = SBMGenerator::fork_generate_sbm(
            unclustered_graph, unclustered_clustering_file, unclustered_sbm_dir, "unclustered", verbose);
            
        if (clustered_pid < 0 || unclustered_pid < 0) {
            std::cerr << "Error forking processes for SBM generation" << std::endl;
        } else {
            // Wait for both processes to complete and load the resulting graphs
            if (verbose) {
                std::cout << "Waiting for SBM generation processes to complete..." << std::endl;
            }
            
            // Only load the clustered SBM into memory
            CSRGraph clustered_sbm = SBMGenerator::wait_and_load_sbm(
                clustered_pid, clustered_sbm_dir, num_threads, verbose);
                
            // For the unclustered SBM, just wait for it to complete
            int status;
            waitpid(unclustered_pid, &status, 0);
            
            if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                if (verbose) {
                    std::cout << "Unclustered SBM generation completed successfully" << std::endl;
                }
                
                // Enforce minimum cluster degrees in the clustered SBM
                if (verbose) {
                    std::cout << "Enforcing minimum cluster degrees in clustered SBM..." << std::endl;
                }

                // Set a minimum degree floor of 1 - clusters will have at least this many internal connections
                uint32_t min_degree_floor = 1;

                // Create a modified version of the clustered SBM with enforced minimum degrees
                CSRGraph modified_clustered_sbm = DegreeEnforcer::enforce_min_cluster_degrees(
                    clustered_sbm, clustered_graph, clustering, clustered_node_map, min_degree_floor, verbose);

                // Create paths
                std::string clustered_sbm_filename = base_name + "_clustered_sbm.tsv";
                std::string unclustered_sbm_filename = unclustered_sbm_dir + "/syn_sbm.tsv";
                std::string merged_sbm_filename = base_name + "_merged_sbm.tsv";

                // Save the modified clustered SBM
                save_graph_edgelist(clustered_sbm_filename, modified_clustered_sbm, verbose);
                
                if (verbose) {
                    std::cout << "Saved clustered SBM to: " << clustered_sbm_filename << std::endl;
                }
                
                // Create the merged SBM by appending the clustered SBM to the unclustered SBM
                if (verbose) {
                    std::cout << "Creating merged SBM by appending clustered SBM to unclustered SBM..." << std::endl;
                }
                
                if (GraphMerger::merge_clustered_to_unclustered(
                        clustered_sbm, unclustered_sbm_filename, merged_sbm_filename, verbose)) {
                    if (verbose) {
                        std::cout << "Successfully created merged SBM: " << merged_sbm_filename << std::endl;
                    }
                } else {
                    std::cerr << "Failed to create merged SBM" << std::endl;
                }
            } else {
                std::cerr << "Unclustered SBM generation failed" << std::endl;
            }
        }
        
        auto sbm_end_time = std::chrono::steady_clock::now();
        auto sbm_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            sbm_end_time - sbm_start_time);
            
        if (verbose) {
            std::cout << "SBM generation time: " << sbm_duration.count() / 1000.0 << " seconds" << std::endl;
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
