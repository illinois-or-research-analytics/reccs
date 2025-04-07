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
};