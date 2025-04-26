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
        // Read file into memory in one operation for faster processing
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
        
        // Process the buffer in parallel to extract unique vertex IDs
        std::vector<std::string> lines;
        std::istringstream stream(std::string(buffer.data(), size));
        std::string line;
        while (std::getline(stream, line)) {
            lines.push_back(line);
        }

        // Use thread-safe set for collecting unique IDs
        std::mutex set_mutex;
        std::set<int> unique_ids;
        
        #pragma omp parallel
        {
            std::set<int> local_ids;
            
            #pragma omp for nowait
            for (size_t i = 0; i < lines.size(); i++) {
                std::istringstream iss(lines[i]);
                int source, target;
                if (iss >> source >> target) {
                    local_ids.insert(source);
                    local_ids.insert(target);
                }
            }
            
            // Merge local sets into the global set
            std::lock_guard<std::mutex> lock(set_mutex);
            unique_ids.insert(local_ids.begin(), local_ids.end());
        }
        
        // Create mapping from original IDs to consecutive indices
        int index = 0;
        for (int id : unique_ids) {
            id_to_index[id] = index++;
        }
        
        // Create a temporary edge list file with remapped IDs
        std::string temp_file = filepath + ".temp";
        std::ofstream temp_out(temp_file);
        
        // Process and write in parallel using local buffers
        #pragma omp parallel
        {
            std::ostringstream local_buffer;
            
            #pragma omp for nowait
            for (size_t i = 0; i < lines.size(); i++) {
                std::istringstream iss(lines[i]);
                int source, target;
                if (iss >> source >> target) {
                    local_buffer << id_to_index[source] << "\t" << id_to_index[target] << "\n";
                }
                
                // Avoid excessive memory usage by periodically flushing large buffers
                if (i % 100000 == 0 && local_buffer.tellp() > 0) {
                    #pragma omp critical
                    {
                        temp_out << local_buffer.str();
                    }
                    local_buffer.str("");
                    local_buffer.clear();
                }
            }
            
            // Final flush of remaining data
            if (local_buffer.tellp() > 0) {
                #pragma omp critical
                {
                    temp_out << local_buffer.str();
                }
            }
        }
        
        temp_out.close();
        
        // Create igraph object
        igraph_t graph_primitive;
        
        // Use igraph's built-in file reading functionality
        FILE* file_ptr = fopen(temp_file.c_str(), "r");
        if (!file_ptr) {
            throw std::runtime_error("Could not open temporary file");
        }
        
        // Read the graph from the file with consecutive IDs
        if (igraph_read_graph_edgelist(&graph_primitive, file_ptr, 0, IGRAPH_UNDIRECTED)) {
            fclose(file_ptr);
            std::remove(temp_file.c_str());
            throw std::runtime_error("Failed to parse graph from file: " + filepath);
        }
        fclose(file_ptr);
        std::remove(temp_file.c_str());
        
        // Create a new Graph object
        Graph graph(graph_primitive, id_to_index);
        
        return graph;
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
