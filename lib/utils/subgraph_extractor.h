class SubgraphExtractor {
public:
    /**
     * Extracts a subgraph from a graph based on a clustering.
     * The subgraph contains only nodes that are assigned to clusters.
     * 
     * @param graph The original graph.
     * @param clustering The clustering of the nodes.
     * @return A subgraph containing only nodes that are assigned to clusters.
     */
    static Graph get_clustered_subgraph(const Graph& graph, Clustering& clustering, std::unordered_map<int, int>& id_to_index) {
        // Remove singletons from the clustering
        clustering.remove_singletons();
        
        // Get the remaining nodes
        std::vector<int> cluster_nodes = clustering.get_nodes();
        
        // Create a new graph for the subgraph
        Graph subgraph = graph.get_induced_subgraph(cluster_nodes);

        return subgraph;
    }

    /**
     * Extracts a subgraph from a graph containing singletons and their neighbors.
     * 
     * @param graph The original graph.
     * @param clustering The clustering of the nodes.
     * @param id_to_index Mapping from original node IDs to indices.
     * @return A subgraph containing singleton nodes and their neighbors.
     */
    static Graph get_outlier_subgraph(
        const Graph& graph,
        Clustering& clustering,
        std::unordered_map<int, int>& id_to_index,
        Clustering& clustering_outlier) {
        // Get the singleton nodes
        std::vector<int> singletons = clustering.get_singletons();
        
        // Get the union of singletons and their neighbors
        std::unordered_set<int> vertices_set(singletons.begin(), singletons.end());
        
        // Add neighbors of singletons to the set
        for (int node : singletons) {
            const auto& neighbors = graph.get_neighbors(node);
            vertices_set.insert(neighbors.begin(), neighbors.end());
        }
        
        // Convert set back to vector
        std::vector<int> vertices(vertices_set.begin(), vertices_set.end());
        
        // Create a new graph for the outlier subgraph with singletons and their neighbors
        Graph outlier_subgraph = graph.get_induced_subgraph(vertices);
        
        // Initialize the clustering_outlier object
        clustering_outlier.clear();
        
        // Add all singletons to a single cluster (cluster -1)
        for (int node : singletons) {
            clustering_outlier.add_node_to_cluster(node, -1);
        }
        
        // For neighbors, preserve their original cluster assignment if they have one
        for (int node : vertices) {
            // Skip if it's a singleton (already added to cluster 0)
            if (std::find(singletons.begin(), singletons.end(), node) != singletons.end()) {
                continue;
            }
            
            // Check if node has a cluster assignment in the original clustering
            int cluster = clustering.get_cluster(node);
            if (cluster != -1) {  // -1 typically indicates no cluster
                // Add to the same cluster in the outlier clustering
                clustering_outlier.add_node_to_cluster(node, cluster);
            }
        }
        
        return outlier_subgraph;
    }
};