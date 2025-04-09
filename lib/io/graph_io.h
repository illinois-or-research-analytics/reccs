#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <set> 
#include <map>

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
    static reccs::Graph load_graph_from_tsv(const std::string& filepath) {
     */
    static Graph load_graph_from_tsv(const std::string& filepath, std::unordered_map<int, int>& id_to_index) {
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
        igraph_t graph_primitive;
        igraph_create(&graph_primitive, &edges, unique_nodes.size(), IGRAPH_UNDIRECTED);

        // Create a new Graph object
        Graph graph(graph_primitive, id_to_index);

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
