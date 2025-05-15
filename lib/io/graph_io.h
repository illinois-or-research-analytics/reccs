#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <omp.h>
#include "../data_structures/graph.h"

/**
 * @brief Utility class for reading and writing graph files
 */
class GraphIO {
public:
    /**
     * @brief Read a graph from a TSV edgelist file
     * 
     * @param filepath Path to the TSV file
     * @param verbose Whether to print verbose information
     * @return Graph The loaded graph
     */
    static Graph read_tsv(const std::string& filepath, bool verbose = false) {
        // Get start time for performance measurement
        auto start_time = std::chrono::high_resolution_clock::now();

        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filepath << std::endl;
            return Graph();
        }
        
        if (verbose) {
            std::cout << "Reading graph from " << filepath << std::endl;
        }
        
        std::vector<int> src, dst;
        std::string line;
        
        // First pass: count the number of lines
        size_t line_count = 0;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;  // Skip empty lines and comments
            line_count++;
        }
        
        if (verbose) {
            std::cout << "Found " << line_count << " edge records" << std::endl;
        }
        
        // Resize vectors to hold all edges
        src.reserve(line_count);
        dst.reserve(line_count);
        
        // Reset file to beginning
        file.clear();
        file.seekg(0, std::ios::beg);
        
        // Read all edges and collect unique vertices
        std::unordered_set<int> unique_vertices;
        
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int source, target;
            
            if (iss >> source >> target) {
                src.push_back(source);
                dst.push_back(target);
                
                unique_vertices.insert(source);
                unique_vertices.insert(target);
            }
        }
        
        file.close();
        
        // Create node mapping from original IDs to continuous indices
        std::unordered_map<int, int> node_map;
        std::vector<int> rev_map;
        rev_map.reserve(unique_vertices.size());
        
        int idx = 0;
        for (int vertex : unique_vertices) {
            node_map[vertex] = idx;
            rev_map.push_back(vertex);
            idx++;
        }
        
        // Create and build the graph with node mapping
        size_t num_edges = src.size();
        
        if (verbose) {
            std::cout << "Creating undirected graph with " << unique_vertices.size() << " vertices and " 
                      << num_edges << " edges" << std::endl;
            
            // Check if IDs are non-continuous
            bool continuous = true;
            for (int i = 0; i < static_cast<int>(unique_vertices.size()); i++) {
                if (unique_vertices.find(i) == unique_vertices.end()) {
                    continuous = false;
                    break;
                }
            }
            
            if (!continuous) {
                std::cout << "Detected non-continuous vertex IDs, creating node mapping" << std::endl;
            }
        }
        
        Graph graph(unique_vertices.size(), num_edges);
        graph.build_from_edges(src, dst, node_map);
        
        if (verbose) {
            std::cout << "Graph loaded successfully." << std::endl;
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
            std::cout << "Loading time: " << duration << " seconds" << std::endl;
            graph.print_stats();
        }
        
        return graph;
    }
    
    /**
     * @brief Write a graph to a TSV edgelist file
     * 
     * @param graph The graph to write
     * @param filepath Path to the output TSV file
     * @param use_original_ids Whether to use original IDs in the output
     * @return bool Success status
     */
    static bool write_tsv(const Graph& graph, const std::string& filepath, bool use_original_ids = true) {
        std::ofstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filepath << " for writing" << std::endl;
            return false;
        }
        
        const auto& src = graph.src();
        const auto& dst = graph.dst();
        
        if (use_original_ids && !graph.node_mapping().empty()) {
            const auto& rev_map = graph.reverse_mapping();
            
            for (size_t i = 0; i < graph.num_edges(); i++) {
                file << rev_map[src[i]] << "\t" << rev_map[dst[i]] << "\n";
            }
        } else {
            for (size_t i = 0; i < graph.num_edges(); i++) {
                file << src[i] << "\t" << dst[i] << "\n";
            }
        }
        
        file.close();
        return true;
    }
    
    /**
     * @brief Read a clustering file (vertex_id cluster_id format)
     * 
     * @param filepath Path to the clustering file
     * @param verbose Whether to print verbose information
     * @return std::vector<int> Cluster assignments for each vertex
     */
    static std::vector<int> read_clustering(const std::string& filepath, bool verbose = false) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open clustering file " << filepath << std::endl;
            return {};
        }
        
        if (verbose) {
            std::cout << "Reading clustering from " << filepath << std::endl;
        }
        
        std::vector<std::pair<int, int>> node_cluster_pairs;
        std::string line;
        
        // Read all vertex-cluster pairs
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int node_id, cluster_id;
            
            if (iss >> node_id >> cluster_id) {
                node_cluster_pairs.emplace_back(node_id, cluster_id);
            }
        }
        
        file.close();
        
        // Find the maximum vertex ID to determine the size of the result array
        int max_vertex_id = -1;
        for (const auto& pair : node_cluster_pairs) {
            max_vertex_id = std::max(max_vertex_id, pair.first);
        }
        
        // Create and fill the clustering array
        std::vector<int> clustering(max_vertex_id + 1, -1);  // -1 indicates no cluster assignment
        
        for (const auto& pair : node_cluster_pairs) {
            clustering[pair.first] = pair.second;
        }
        
        if (verbose) {
            std::cout << "Clustering loaded successfully for " << node_cluster_pairs.size() << " vertices" << std::endl;
            
            // Count the number of unique clusters
            std::unordered_set<int> unique_clusters;
            for (const auto& pair : node_cluster_pairs) {
                unique_clusters.insert(pair.second);
            }
            std::cout << "Found " << unique_clusters.size() << " unique clusters" << std::endl;
        }
        
        return clustering;
    }
    
    /**
     * @brief Read a clustering file and map to internal vertex IDs
     * 
     * @param filepath Path to the clustering file
     * @param graph The graph with node mapping
     * @param verbose Whether to print verbose information
     * @return std::vector<int> Cluster assignments for each vertex (using internal IDs)
     */
    static std::vector<int> read_clustering_mapped(const std::string& filepath, 
                                                 const Graph& graph, 
                                                 bool verbose = false) {
        std::vector<std::pair<int, int>> node_cluster_pairs;
        std::ifstream file(filepath);
        
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open clustering file " << filepath << std::endl;
            return {};
        }
        
        if (verbose) {
            std::cout << "Reading clustering from " << filepath << std::endl;
        }
        
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream iss(line);
            int node_id, cluster_id;
            
            if (iss >> node_id >> cluster_id) {
                node_cluster_pairs.emplace_back(node_id, cluster_id);
            }
        }
        
        file.close();
        
        // Create clustering array for internal vertex IDs
        std::vector<int> clustering(graph.num_vertices(), -1);
        const auto& node_map = graph.node_mapping();
        
        for (const auto& pair : node_cluster_pairs) {
            int original_id = pair.first;
            int cluster_id = pair.second;
            
            auto it = node_map.find(original_id);
            if (it != node_map.end()) {
                int internal_id = it->second;
                clustering[internal_id] = cluster_id;
            }
        }
        
        if (verbose) {
            // Count vertices with cluster assignments
            int assigned = 0;
            for (int c : clustering) {
                if (c != -1) assigned++;
            }
            
            std::cout << "Mapped clustering to " << assigned << " vertices (out of " 
                     << graph.num_vertices() << ")" << std::endl;
            
            // Count unique clusters
            std::unordered_set<int> unique_clusters;
            for (int c : clustering) {
                if (c != -1) unique_clusters.insert(c);
            }
            std::cout << "Found " << unique_clusters.size() << " unique clusters" << std::endl;
        }
        
        return clustering;
    }
    
    /**
     * @brief Write a clustering to a file (vertex_id cluster_id format)
     * 
     * @param clustering The clustering assignments (using internal IDs)
     * @param graph The graph with node mapping
     * @param filepath Path to the output file
     * @param use_original_ids Whether to use original IDs in the output
     * @return bool Success status
     */
    static bool write_clustering(const std::vector<int>& clustering, 
                               const Graph& graph,
                               const std::string& filepath,
                               bool use_original_ids = true) {
        std::ofstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Error: Unable to open file " << filepath << " for writing" << std::endl;
            return false;
        }
        
        if (use_original_ids && !graph.node_mapping().empty()) {
            const auto& rev_map = graph.reverse_mapping();
            
            for (size_t i = 0; i < clustering.size(); i++) {
                if (clustering[i] != -1) {  // Only write assigned vertices
                    file << rev_map[i] << "\t" << clustering[i] << "\n";
                }
            }
        } else {
            for (size_t i = 0; i < clustering.size(); i++) {
                if (clustering[i] != -1) {  // Only write assigned vertices
                    file << i << "\t" << clustering[i] << "\n";
                }
            }
        }
        
        file.close();
        return true;
    }
};
