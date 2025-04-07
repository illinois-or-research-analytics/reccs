#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <set> 
#include <map>

#include "data_structures/graph.h"
#include "data_structures/clustering.h"

class graph_io {
public:
    /**
     * Loads a graph from a TSV file containing an edgelist.
     * Each line contains: source_id \t target_id
     * 
     * @param filepath Path to the TSV file
    static reccs::Graph load_graph_from_tsv(const std::string& filepath) {
     */
    static igraph_t load_graph_from_tsv(const std::string& filepath, std::map<int, int>& id_to_index) {
        // Vectors to store edges
        std::vector<int> from_nodes;
        std::vector<int> to_nodes;
        std::set<int> unique_nodes;  // To track all unique node IDs
        
        // Read the TSV file
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filepath << std::endl;
            exit(EXIT_FAILURE);
        }
        
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int from, to;
            
            // Parse TSV line (tab-separated values)
            if (!(iss >> from >> to)) {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
            
            // Store edge endpoints
            from_nodes.push_back(from);
            to_nodes.push_back(to);
            
            // Keep track of unique node IDs
            unique_nodes.insert(from);
            unique_nodes.insert(to);
        }
        file.close();
        
        // Create mapping from original IDs to contiguous indices
        int index = 0;
        for (int id : unique_nodes) {
            id_to_index[id] = index++;
        }
        
        // Create igraph edge vector with mapped indices
        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, from_nodes.size() * 2);
        
        for (size_t i = 0; i < from_nodes.size(); i++) {
            VECTOR(edges)[2*i] = id_to_index[from_nodes[i]];
            VECTOR(edges)[2*i+1] = id_to_index[to_nodes[i]];
        }
        
        // Create the graph
        igraph_t graph;
        igraph_create(&graph, &edges, unique_nodes.size(), IGRAPH_UNDIRECTED);
        
        // Clean up
        igraph_vector_int_destroy(&edges);

        return graph;
    }
    
    /**
     * Loads a clustering from a TSV file.
     * Each line contains: node_id \t cluster_id
     * 
     * @param filepath Path to the TSV file
    static reccs::Clustering load_clustering_from_tsv(const std::string& filepath) {
     */
    static Clustering load_clustering_from_tsv(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }
        
        Clustering clustering;
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int node_id;
            int cluster_id;
            if (!(iss >> node_id >> cluster_id)) {
                continue; // Skip malformed lines
            }
            clustering.add_node_to_cluster(node_id, cluster_id);
        }

        return clustering;
    }
};
