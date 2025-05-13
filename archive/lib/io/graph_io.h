#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <set> 
#include <map>

// Include OpenMP for parallel processing
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <mutex>

#include "data_structures/graph.h"
#include "data_structures/clustering.h"
#include "data_structures/edge_iterator.h"

class graph_io {
public:
    /**
     * Loads a graph from a TSV file containing an edgelist.
     * Each line contains: source_id \t target_id
     * 
     * @param filepath Path to the TSV file
     */
    static Graph load_graph_from_tsv(const std::string& filepath, 
                                      std::unordered_map<int, int>& id_to_index) {
        // Initialize igraph
        igraph_t graph;
        
        // Cast filepath to a C string
        const char* c_filepath = filepath.c_str();

        // Open the input file
        FILE* instream = fopen(c_filepath, "r");
        if (!instream) {
            std::cerr << "Failed to open file: " << filepath << std::endl;
            exit(1);
        }
        
        // No predefined vertex names
        igraph_strvector_t predefnames;
        igraph_strvector_init(&predefnames, 0);
        
        // Read the graph using the NCOL format
        // Parameters:
        // - graph: output graph
        // - instream: input file
        // - predefnames: predefined vertex names (empty in this case)
        // - names: true to use symbolic names in the file
        // - weights: IGRAPH_ADD_WEIGHTS_IF_PRESENT to add weights if present
        // - directed: false for undirected graph
        igraph_error_t err = igraph_read_graph_ncol(
            &graph, 
            instream, 
            &predefnames,
            true,  // Use symbolic names from the file
            IGRAPH_ADD_WEIGHTS_IF_PRESENT,  // Add weights if present
            false  // Undirected graph
        );
        
        // Close the file
        fclose(instream);
        
        if (err != IGRAPH_SUCCESS) {
            std::cerr << "Failed to read graph. Error code: " << err << std::endl;
            igraph_strvector_destroy(&predefnames);
            exit(1);
        }

        // Destroy the predefined names vector
        // igraph_strvector_destroy(&predefnames);

        // Convert igraph to our Graph class
        Graph g(graph);
    }
    
    /**
     * Loads a clustering from a TSV file.
     * Each line contains: node_id \t cluster_id
     * 
     * @param filepath Path to the TSV file
     */
    static Clustering load_clustering_from_tsv(const std::string& filepath) {
        // Read file into memory in one operation
        std::ifstream file(filepath, std::ios::ate);  // Open at end to get file size
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }
        
        // Pre-allocate memory for the file content
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        std::vector<char> buffer(size);
        file.read(buffer.data(), size);
        file.close();
        
        // Process the buffer
        Clustering clustering;
        std::unordered_map<int, std::vector<int>> tmp_clusters;
        
        const char* p = buffer.data();
        const char* end = p + size;
        
        // Count number of lines to pre-allocate vectors
        size_t line_count = std::count(p, end, '\n');
        
        while (p < end) {
            // Parse node_id
            int node_id = 0;
            while (p < end && *p >= '0' && *p <= '9') {
                node_id = node_id * 10 + (*p - '0');
                p++;
            }
            
            // Skip whitespace
            while (p < end && (*p == ' ' || *p == '\t')) p++;
            
            // Parse cluster_id
            int cluster_id = 0;
            while (p < end && *p >= '0' && *p <= '9') {
                cluster_id = cluster_id * 10 + (*p - '0');
                p++;
            }
            
            // Add to clustering
            if (!tmp_clusters.count(cluster_id)) {
                tmp_clusters[cluster_id].reserve(line_count / 10);  // Rough estimate
            }
            tmp_clusters[cluster_id].push_back(node_id);
            
            // Skip to next line
            while (p < end && *p != '\n') p++;
            if (p < end) p++;  // Skip newline
        }
        
        // Build the clustering from temp structure
        for (const auto& [cluster_id, nodes] : tmp_clusters) {
            for (int node_id : nodes) {
                clustering.add_node_to_cluster(node_id, cluster_id);
            }
        }
        
        return clustering;
    }

    /**
     * Writes a graph to a TSV file containing an edgelist.
     * Each line contains: source_id \t target_id
     * 
     * @param graph The graph to write
     * @param output_file Path to the output TSV file
     */
    static void write_graph_to_tsv(const Graph& graph, const std::string& output_file) {
        // Open the output file
        EdgeIterator edge_iterator(graph);
        std::ofstream output(output_file);
        if (!output.is_open()) {
            std::cerr << "Error opening output file: " << output_file << std::endl;
            exit(1);
        }
        
        // Write the edges to the file
        for (edge_iterator.reset(); edge_iterator.has_next(); edge_iterator.next()) {
            igraph_integer_t from, to;
            edge_iterator.get(from, to);
            
            output << from << "\t" << to << "\n";
        }

        // Close the output file
        output.close();
    }

    /**
     * Writes a clustering to a TSV file.
     * Each line contains: node_id \t cluster_id
     * 
     * @param clustering The clustering to write
     * @param output_file Path to the output TSV file
     */
    static void write_clustering_to_tsv(const Clustering& clustering, const std::string& output_file) {
        // Open the output file
        std::ofstream output(output_file);
        if (!output.is_open()) {
            std::cerr << "Error opening output file: " << output_file << std::endl;
            exit(1);
        }
        
        // Write the clustering to the file
        std::unordered_map<int, std::vector<int>> clusters = clustering.get_clusters();
        for (const auto& cluster : clusters) {
            for (int node_id : cluster.second) {
                output << node_id << "\t" << cluster.first << "\n";
            }
        }

        // Close the output file
        output.close();
    }
};
