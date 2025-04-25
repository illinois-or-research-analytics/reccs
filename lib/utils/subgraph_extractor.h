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
    static Graph get_outlier_subgraph(const Graph& graph, Clustering& clustering, std::unordered_map<int, int>& id_to_index) {
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
        
        return outlier_subgraph;
    }
};