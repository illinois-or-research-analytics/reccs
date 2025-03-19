#pragma once

#include <iostream>
#include "data_structures/graph.h"
#include "data_structures/clustering.h"

// namespace reccs {

class graph_io {
public:
    /**
     * Loads a graph from a TSV file containing an edgelist.
     * Each line contains: source_id \t target_id
     * 
     * @param filepath Path to the TSV file
     * @return Graph object constructed from the edgelist
     */
    static Graph load_graph_from_tsv(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open file: " + filepath);
        }
        
        Graph graph;
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            node_id_t src, dst;
            if (!(iss >> src >> dst)) {
                continue; // Skip malformed lines
            }
            graph.add_edge(src, dst);
        }
        
        return graph;
    }
    
    /**
     * Loads a clustering from a TSV file.
     * Each line contains: node_id \t cluster_id
     * 
     * @param filepath Path to the TSV file
     * @return Clustering object mapping nodes to clusters
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
            node_id_t node_id;
            cluster_id_t cluster_id;
            if (!(iss >> node_id >> cluster_id)) {
                continue; // Skip malformed lines
            }
            clustering.add_node_to_cluster(node_id, cluster_id);
        }
        
        return clustering;
    }
};

// } // namespace reccs
