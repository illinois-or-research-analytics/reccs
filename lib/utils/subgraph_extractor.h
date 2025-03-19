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
    static Graph<int> extract_clustered_subgraph(const Graph<int>& graph, Clustering& clustering) {
        clustering.remove_singletons();

        Graph<int> subgraph = graph.induced_subgraph(clustering.get_nodes());
        return subgraph;
    }
};