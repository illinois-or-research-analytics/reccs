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
    static igraph_t get_clustered_subgraph(const igraph_t& graph, Clustering& clustering, std::map<int, int>& id_to_index) {
        // Remove singletons from the clustering
        clustering.remove_singletons();
        
        // Get the remaining nodes
        std::vector<int> cluster_nodes = clustering.get_nodes();
        
        // Create a vertex selector from the nodes vector
        igraph_vs_t vs;
        igraph_vector_int_t nodes_vec;
        igraph_vector_int_init(&nodes_vec, cluster_nodes.size());
        
        // Fill the igraph vector with node IDs from the vector
        for (size_t i = 0; i < cluster_nodes.size(); i++) {
            VECTOR(nodes_vec)[i] = id_to_index[cluster_nodes[i]];
        }
        
        // Create vertex selector from the vector
        igraph_vs_vector(&vs, &nodes_vec);
        
        // Get the induced subgraph
        igraph_t subgraph;
        igraph_induced_subgraph(&graph, &subgraph, vs, IGRAPH_SUBGRAPH_AUTO);
        
        // Update id_to_index map with the mapping from original IDs to new indices
        id_to_index.clear();
        for (size_t i = 0; i < cluster_nodes.size(); i++) {
            id_to_index[cluster_nodes[i]] = i;
        }
        
        // Clean up
        igraph_vs_destroy(&vs);
        igraph_vector_int_destroy(&nodes_vec);
        
        return subgraph;
    }
};