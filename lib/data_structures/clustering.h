#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <unordered_map>
#include <vector>
#include <stdexcept>

// namespace reccs {

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
};

// } // namespace reccs

#endif // CLUSTERING_H