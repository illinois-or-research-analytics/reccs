#ifndef DEG_SEQ_MATCHING_V1_H
#define DEG_SEQ_MATCHING_V1_H

#include <iostream>
#include <string>
#include <filesystem>
#include <cstdlib>
#include "../data_structures/graph.h"
#include "../io/g_io.h"

/**
 * Match the degree sequence of a clustered SBM graph using a simple greedy approach.
 */
void match_degree_sequence_v1(
    Graph& g,
    std::string empirical_graph_path,
    std::string clustering_path,
    const std::string& temp_dir,
    std::string output_file) {

    // Create temp directory if it doesn't exist
    std::filesystem::create_directories(temp_dir);
    std::filesystem::create_directories(temp_dir + "/v1_output");

    std::cout << "Starting degree sequence matching for graph " << g.id << std::endl;

    // Step 1: Save current graph state
    std::string current_edges_file = temp_dir + "/current_edges.tsv";
    save_graph_edgelist(current_edges_file, g, true);
    std::cout << "Saved current graph to: " << current_edges_file << std::endl;
    
    // Step 2: Execute Python matching script
    std::string match_command = "python3 extlib/correct_degree.py";
    match_command += " --input-edgelist " + current_edges_file; // Current synthetic graph
    match_command += " --ref-edgelist " + empirical_graph_path; // Target
    match_command += " --ref-clustering " + clustering_path; // Clustering file
    match_command += " --output-folder " + temp_dir + "/v1_output"; // Output directory

    std::cout << "Executing degree sequence matching command: " << match_command << std::endl;

    // Execute the Python script
    int result = std::system(match_command.c_str());

    if (result != 0) {
        std::cerr << "Error: Degree sequence matching script failed with return code " << result << std::endl;
        return;
    }

    // Rename the output file
    std::string new_edges_filename = temp_dir + "/v1_output/degcorr_edge.tsv";
    
    // Append the output file to the existing current edges tsv file
    std::ofstream output_stream(output_file);
    std::ifstream new_edges_stream(new_edges_filename);
    std::ifstream current_edges_stream(current_edges_file);

    // Append current edges first
    output_stream << current_edges_stream.rdbuf();
    current_edges_stream.close();

    // Then append new edges
    output_stream << new_edges_stream.rdbuf();
    new_edges_stream.close();

    output_stream.close();

    std::cout << "Degree sequence matching completed successfully. Output saved to: " << output_file << std::endl;
}

#endif // DEG_SEQ_MATCHING_V1_H
