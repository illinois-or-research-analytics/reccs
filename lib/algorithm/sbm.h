#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <random>

#include "io/graph_io.h"

// New class to compute and work with SBM
class StochasticBlockModel {
private:
    std::vector<std::vector<double>> block_matrix;
    int num_blocks;

public:
    // Compute SBM from a graph and clustering
    StochasticBlockModel(const Graph<int>& graph, const Clustering& clustering) {
        num_blocks = clustering.get_clusters().size();
        
        // Initialize block matrix with zeros
        block_matrix.resize(num_blocks, std::vector<double>(num_blocks, 0.0));
        
        // Count nodes in each block
        std::vector<int> block_sizes(num_blocks, 0);
        for (const auto& cluster_pair : clustering.get_clusters()) {
            int cluster_id = cluster_pair.first;
            block_sizes[cluster_id] = cluster_pair.second.size();
        }
        
        // Count edges between each pair of blocks
        std::vector<std::vector<int>> edge_counts(num_blocks, std::vector<int>(num_blocks, 0));
        
        // Iterate over all nodes using the NodeIterator and count edges between blocks
        for (auto node_it = graph.begin(); node_it != graph.end(); ++node_it) {
            int node_i = *node_it;
            int block_i = clustering.get_cluster(node_i);
            
            for (const auto& neighbor : graph.get_neighbors(node_i)) {
                int block_j = clustering.get_cluster(neighbor);
                edge_counts[block_i][block_j]++;
            }
        }
        
        // Compute probabilities
        for (int i = 0; i < num_blocks; i++) {
            for (int j = 0; j < num_blocks; j++) {
                // Calculate maximum possible edges between blocks
                int max_edges;
                if (i == j) {
                    // Within a block: n(n-1)/2 possible edges
                    max_edges = block_sizes[i] * (block_sizes[i] - 1) / 2;
                } else {
                    // Between blocks: n_i * n_j possible edges
                    max_edges = block_sizes[i] * block_sizes[j];
                }
                
                // Avoid division by zero
                if (max_edges > 0) {
                    // Each edge is counted twice in undirected graphs, so divide by 2
                    block_matrix[i][j] = static_cast<double>(edge_counts[i][j]) / (2 * max_edges);
                }
            }
        }
    }
    
    // Get the block matrix
    const std::vector<std::vector<double>>& get_block_matrix() const {
        return block_matrix;
    }
    
    // Print the block matrix
    void print_block_matrix() const {
        std::cout << "SBM Block Matrix:" << std::endl;
        for (int i = 0; i < num_blocks; i++) {
            for (int j = 0; j < num_blocks; j++) {
                std::cout << block_matrix[i][j] << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    // Generate a new graph from this SBM
    Graph<int> generate_graph(int num_nodes, const std::vector<int>& block_assignments) {
        Graph<int> new_graph;
        
        // Add nodes
        for (int i = 0; i < num_nodes; i++) {
            new_graph.add_node(i);
        }
        
        // Precompute all bernoulli distributions
        std::vector<std::vector<std::bernoulli_distribution>> edge_dists(num_blocks, 
            std::vector<std::bernoulli_distribution>(num_blocks));
        
        for (int i = 0; i < num_blocks; i++) {
            for (int j = 0; j < num_blocks; j++) {
                edge_dists[i][j] = std::bernoulli_distribution(block_matrix[i][j]);
            }
        }
        
        // Add edges according to block probabilities
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for (int i = 0; i < num_nodes; i++) {
            int block_i = block_assignments[i];
            
            for (int j = i + 1; j < num_nodes; j++) {
                int block_j = block_assignments[j];
                
                // Use the precomputed distribution
                if (edge_dists[block_i][block_j](gen)) {
                    new_graph.add_edge(i, j);
                }
            }
        }
        
        return new_graph;
    }
};

