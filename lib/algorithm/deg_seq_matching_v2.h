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
    const std::string& temp_dir,
    std::string output_file) {

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

    // Call the tsv union script to merge results
    std::string append_command = "python3 extlib/tsv_union.py";
    append_command += " " + temp_dir + "/current_edges.tsv"; // Current graph
    append_command += " " + temp_dir + "/sbm_output/new_edges.tsv"; // New edges from SBM
    append_command += " " + output_file; // Final output file

    std::cout << "Executing TSV union command: " << append_command << std::endl;

    // Execute the append command
    result = std::system(append_command.c_str());

    if (result != 0) {
        std::cerr << "Error: TSV union script failed with return code " << result << std::endl;
        return;
    }
}

#endif // DEG_SEQ_MATCHING_V2_H
