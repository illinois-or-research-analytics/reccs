#include <iostream>
#include <chrono>
#include <omp.h>
#include <string>
#include <iomanip>
#include "../lib/data_structures/graph.h"
#include "../lib/data_structures/clustering.h"
#include "../lib/io/graph_io.h"
#include "../lib/io/cluster_io.h"
#include "../lib/algorithms/graph_splitter.h"

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name << " <edgelist.tsv> [options]" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -t <num_threads>  Number of threads to use (default: hardware concurrency)" << std::endl;
    std::cerr << "  -v                Verbose mode: print detailed progress information" << std::endl;
    std::cerr << "  -c <clusters.tsv> Load clusters from TSV file" << std::endl;
    std::cerr << "  -o <output_prefix> Prefix for output files (default: input filename without extension)" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    // Parse command line arguments
    std::string graph_filename;
    std::string cluster_filename;
    std::string output_prefix;
    int num_threads = std::thread::hardware_concurrency();
    bool verbose = false;
    
    // First argument is the graph filename
    graph_filename = argv[1];
    
    // Default output prefix is input filename without extension
    size_t dot_pos = graph_filename.find_last_of('.');
    size_t slash_pos = graph_filename.find_last_of('/');
    if (dot_pos != std::string::npos) {
        if (slash_pos != std::string::npos && slash_pos > dot_pos) {
            // No extension in the filename itself, use whole filename
            output_prefix = graph_filename;
        } else {
            // Remove extension
            output_prefix = graph_filename.substr(0, dot_pos);
        }
    } else {
        // No extension
        output_prefix = graph_filename;
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
        } else if (arg == "-o" && i + 1 < argc) {
            output_prefix = argv[++i];
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
        std::cout << "Output prefix: " << output_prefix << std::endl;
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
                }
            }
        }
        
        // Split the graph (always done when clustering is loaded)
        auto split_start_time = std::chrono::steady_clock::now();
        
        auto [clustered_graph, unclustered_graph] = 
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
        
        // Save the split graphs
        auto save_start_time = std::chrono::steady_clock::now();
        
        std::string clustered_filename = output_prefix + "_clustered.tsv";
        save_graph_edgelist(clustered_filename, clustered_graph, verbose);
        
        std::string unclustered_filename = output_prefix + "_unclustered.tsv";
        save_graph_edgelist(unclustered_filename, unclustered_graph, verbose);
        
        auto save_end_time = std::chrono::steady_clock::now();
        auto save_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            save_end_time - save_start_time);
            
        if (verbose) {
            std::cout << "Graph saving time: " << save_duration.count() / 1000.0 << " seconds" << std::endl;
            std::cout << "Saved clustered graph to: " << clustered_filename << std::endl;
            std::cout << "Saved unclustered graph to: " << unclustered_filename << std::endl;
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
