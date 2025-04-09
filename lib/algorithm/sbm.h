#include <ctime>

/**
 * @brief Stochastic Block Model (SBM) for generating random graphs.
 */
class StochasticBlockModel {
private:
    igraph_matrix_t block_matrix;
    igraph_vector_int_t cluster_sizes;
    int num_blocks;

public:
    /**
     * @brief Constructor for StochasticBlockModel.
     * 
     * @param graph The input graph.
     * @param clustering The clustering of the nodes.
     */
    StochasticBlockModel(const Graph& graph, const Clustering& clustering) {
        // Initialize the block matrix
        num_blocks = clustering.get_num_clusters();
        igraph_matrix_init(&block_matrix, num_blocks, num_blocks);
        igraph_vector_int_init(&cluster_sizes, num_blocks);

        for (int i = 0; i < num_blocks; i++) {
            VECTOR(cluster_sizes)[i] = clustering.get_cluster_size(i);
        }

        // Fill the block matrix
        EdgeIterator edge_it(graph);
        for (edge_it.reset(); edge_it.has_next(); edge_it.next()) {
            igraph_integer_t from, to;
            edge_it.get(from, to);

            int cluster_from = clustering.get_cluster(from);
            int cluster_to = clustering.get_cluster(to);
            
            if (cluster_from == cluster_to) {
                MATRIX(block_matrix, cluster_from, cluster_to) += 1;
            } else {
                MATRIX(block_matrix, cluster_from, cluster_to) += 1;
                MATRIX(block_matrix, cluster_to, cluster_from) += 1;
            }
        }

        // Normalize the block matrix
        for (int i = 0; i < num_blocks; i++) {
            for (int j = 0; j < num_blocks; j++) {
                if (VECTOR(cluster_sizes)[i] > 0 && VECTOR(cluster_sizes)[j] > 0) {
                    if (i == j) {
                        // Within-cluster: divide by the number of possible edges
                        int n = VECTOR(cluster_sizes)[i];
                        MATRIX(block_matrix, i, j) /= (n * (n-1) / 2.0);
                    } else {
                        // Between-cluster: divide by the product of sizes
                        MATRIX(block_matrix, i, j) /= (VECTOR(cluster_sizes)[i] * VECTOR(cluster_sizes)[j]);
                    }
                }
            }
        }
    }
    
    /**
     * @brief Get the block matrix.
     */
    const igraph_matrix_t* get_block_matrix() const {
        return &block_matrix;
    }
    
    /**
     * @brief Get the number of blocks in the SBM.
     * 
     * @return The number of blocks.
     */
    void print_block_matrix() const {
        std::cout << "SBM Block Matrix:" << std::endl;
        for (int i = 0; i < num_blocks; i++) {
            for (int j = 0; j < num_blocks; j++) {
                std::cout << MATRIX(block_matrix, i, j) << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    /**
     * @brief Generate a new graph from the SBM.
     * 
     * @return A new graph generated from the SBM.
     */
    Graph generate_graph() {
        // Create an igraph_t object
        igraph_t g;
        
        // Seed the random number generator with current time
        igraph_rng_seed(igraph_rng_default(), time(NULL));
                
        // Generate a graph according to the SBM model
        igraph_error_t err = igraph_sbm_game(&g, 
                                            igraph_vector_int_sum(&cluster_sizes), 
                                            &block_matrix, 
                                            &cluster_sizes, 
                                            /* directed = */ false, 
                                            /* loops = */ false);
                                            
        if (err != IGRAPH_SUCCESS) {
            throw std::runtime_error("Failed to generate SBM graph");
        }

        // Create a mapping from original IDs to contiguous indices. Make it an identity mapping
        std::unordered_map<int, int> id_to_index;
        for (int i = 0; i < igraph_vcount(&g); i++) {
            id_to_index[i] = i;
        }

        // Convert the igraph_t to our Graph object
        Graph new_graph(g, id_to_index);

        // No need to remove loops, but remove multiple edges
        new_graph.cleanup(true, false); // Loops are already removed in the igraph_sbm_game function

        return new_graph;
    }
};
