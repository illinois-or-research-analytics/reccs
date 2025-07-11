#ifndef DEG_SEQ_MATCHING_V2_H
#define DEG_SEQ_MATCHING_V2_H

#include <vector>
#include <memory>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdlib>
#include "../data_structures/graph.h"
#include "../data_structures/node_degree.h"
#include "../io/g_io.h"

/**
 * Matches the degree sequence of a graph using hybrid SBM generation.
 * Uses the original working strategy: compares current graph to empirical graph.
 */
void match_degree_sequence_v2(
    Graph& g, 
    std::string empirical_graph_path, 
    std::string clustering_path, 
    const std::string& temp_dir) {

    // Create temp directory if it doesn't exist
    std::filesystem::create_directories(temp_dir);
    std::filesystem::create_directories(temp_dir + "/sbm_output");

    std::cout << "Starting hybrid SBM-based degree sequence matching for graph " << g.id << std::endl;

    // Step 1: Save current graph state
    std::string current_edges_file = temp_dir + "/current_edges.tsv";
    save_graph_edgelist(current_edges_file, g, true);
    std::cout << "Saved current graph to: " << current_edges_file << std::endl;

    // Step 2: Execute Python hybrid SBM script
    std::string sbm_command = "python3 extlib/gen_SBM_degseq.py";
    sbm_command += " -f " + current_edges_file;        // Current synthetic graph
    sbm_command += " -ef " + empirical_graph_path;     // Target empirical graph
    sbm_command += " -c " + clustering_path;           // Clustering file
    sbm_command += " -o " + temp_dir + "/sbm_output";  // Output directory
    sbm_command += " -v";                              // Verbose output
    
    std::cout << "Executing hybrid SBM command: " << sbm_command << std::endl;
    
    // Execute the Python script
    int result = std::system(sbm_command.c_str());
    
    if (result != 0) {
        std::cerr << "Error: Hybrid SBM script failed with return code " << result << std::endl;
        return;
    }

    std::cout << "Hybrid SBM script completed successfully" << std::endl;

    // Step 3: Load generated edges and add them to the graph
    std::string new_edges_file = temp_dir + "/sbm_output/new_edges.tsv";
    
    if (!std::filesystem::exists(new_edges_file)) {
        std::cerr << "Warning: No new edges file generated: " << new_edges_file << std::endl;
        return;
    }

    // Load and add edges using existing infrastructure
    std::ifstream file(new_edges_file);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open edges file " << new_edges_file << std::endl;
        return;
    }

    std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string source_str, target_str;
        
        if (std::getline(iss, source_str, '\t') && std::getline(iss, target_str)) {
            try {
                uint32_t source = std::stoull(source_str);
                uint32_t target = std::stoull(target_str);
                edges_to_add.emplace_back(source, target);
            } catch (const std::exception& e) {
                std::cerr << "Warning: Could not parse edge: " << line << " - " << e.what() << std::endl;
            }
        }
    }
    
    file.close();
    
    // Add edges using batch method
    size_t initial_edges = g.num_edges;
    add_edges_batch(g, edges_to_add);
    uint32_t edges_added = g.num_edges - initial_edges;
    
    std::cout << "[Graph " << g.id << "]: Successfully added " << edges_added 
              << " edges using hybrid SBM generation" << std::endl;
}

#endif // DEG_SEQ_MATCHING_V2_H
