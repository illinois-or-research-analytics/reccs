#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <unordered_map>
#include <vector>
#include <stdexcept>

/**
 * @brief Class that efficiently stores mappings from node IDs to cluster IDs.
 */
class Clustering {
private:
    std::unordered_map<int, int> node_to_cluster;  // Maps node_id to cluster_id

public:
    /**
     * @brief Default constructor.
     */
    Clustering() = default;

    /**
     * @brief Get all nodes that have been assigned to clusters.
     * 
     * @return A vector containing all node IDs with cluster assignments.
     */
    std::vector<int> get_nodes() const {
        std::vector<int> nodes;
        nodes.reserve(node_to_cluster.size());
        for (const auto& [node_id, _] : node_to_cluster) {
            nodes.push_back(node_id);
        }
        return nodes;
    }

    /**
     * @brief Adds a node to a specific cluster.
     * 
     * @param node_id The ID of the node to add.
     * @param cluster_id The ID of the cluster to add the node to.
     */
    void add_node_to_cluster(int node_id, int cluster_id) {
        node_to_cluster[node_id] = cluster_id;
    }

    /**
     * @brief Get the cluster ID for a given node.
     * 
     * @param node_id The ID of the node to query.
     * @return The cluster ID of the node.
     * @throws std::out_of_range if the node is not assigned to any cluster.
     */
    int get_cluster(int node_id) const {
        return node_to_cluster.at(node_id);
    }

    /**
     * @brief Check if a node has been assigned to a cluster.
     * 
     * @param node_id The ID of the node to check.
     * @return true if the node is assigned to a cluster, false otherwise.
     */
    bool has_node(int node_id) const {
        return node_to_cluster.find(node_id) != node_to_cluster.end();
    }

    /**
     * @brief Get the number of nodes assigned to clusters.
     * 
     * @return The number of nodes with cluster assignments.
     */
    size_t size() const {
        return node_to_cluster.size();
    }

    /**
     * @brief Clear all node-to-cluster mappings.
     */
    void clear() {
        node_to_cluster.clear();
    }

    /**
     * @brief Get the number of unique clusters.
     * 
     * @return The number of unique clusters.
     */
    int get_num_clusters() {
        std::unordered_set<int> clusters;
        for (const auto& [node, cluster] : node_to_cluster) {
            clusters.insert(cluster);
        }
        return clusters.size();
    }

    /**
     * @brief Get all singleton nodes (nodes in their own single-node cluster).
     * 
     * @return A vector of node IDs that are singletons.
     */
    std::vector<int> get_singletons() const {
        std::vector<int> singletons;
        std::unordered_map<int, std::vector<int>> cluster_to_nodes;
        
        // Group nodes by cluster
        for (const auto& [node_id, cluster_id] : node_to_cluster) {
            cluster_to_nodes[cluster_id].push_back(node_id);
        }
        
        // Find clusters with only one node
        for (const auto& [cluster_id, nodes] : cluster_to_nodes) {
            if (nodes.size() == 1) {
                singletons.push_back(nodes[0]);
            }
        }
        
        return singletons;
    }

    /**
     * @brief Remove all singleton nodes from the clustering.
     */
    void remove_singletons() {
        for (int node_id : get_singletons()) {
            node_to_cluster.erase(node_id);
        }
    }

    /**
     * @brief Get all clusters and their member nodes.
     * 
     * @return An unordered map mapping cluster IDs to vectors of node IDs.
     */
    std::unordered_map<int, std::vector<int>> get_clusters() const {
        std::unordered_map<int, std::vector<int>> cluster_to_nodes;
        
        // Group nodes by cluster
        for (const auto& [node_id, cluster_id] : node_to_cluster) {
            cluster_to_nodes[cluster_id].push_back(node_id);
        }
        
        return cluster_to_nodes;
    }
    
    /**
     * @brief Return a vector for cluster assignments
     * This is for SBM input
     * 
     * @return Cluster assignments vector
     */
    std::vector<int> get_block_assignments() const {
        std::vector<int> assignments;
        assignments.reserve(node_to_cluster.size());
        
        for (const auto& [node_id, cluster_id] : node_to_cluster) {
            assignments.emplace_back(cluster_id);
        }
        
        return assignments;
    }

    /**
     * @brief Get the subgraph corresponding to a specific cluster.
     * 
     * @param graph The original graph.
     * @param cluster_id The ID of the cluster to extract.
     * @return A subgraph containing only nodes in the specified cluster.
     */
    Graph<int> get_subgraph(const Graph<int>& graph, int cluster_id) const {
        std::vector<int> cluster_nodes = get_clusters()[cluster_id];
        Graph<int> subgraph = graph.induced_subgraph(cluster_nodes);
        return subgraph;
    }
};

#endif // CLUSTERING_H