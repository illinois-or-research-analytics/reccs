#include "data_structures/graph.h"

/**
* @brief Class to iterate through edges in the graph.
*/
class EdgeIterator {
private:
    const Graph& graph;  // Reference to the graph
    const igraph_t& igraph;  // Reference to the igraph object
    igraph_es_t es;
    igraph_eit_t eit;
    
public:
    /**
     * @brief Constructor for EdgeIterator.
     * 
     * @param graph The graph to iterate over.
     */
    EdgeIterator(const Graph& graph) : 
        graph(graph), 
        igraph(graph.get_igraph()) // Initialize the reference in the initialization list
    {
        // Initialize edge selector and iterator
        igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
        igraph_eit_create(&igraph, es, &eit);

        next();
    }

    /**
     * @brief Destructor for EdgeIterator.
     */
    ~EdgeIterator() {
        igraph_eit_destroy(&eit);
        igraph_es_destroy(&es);
    }

    /**
     * @brief Reset the iterator to the beginning.
     */
    void reset() {
        IGRAPH_EIT_RESET(eit);
    }

    /**
     * @brief Check if the iterator has reached the end.
     * 
     * @return bool True if there are more edges, false if at the end.
     */
    bool has_next() const {
        return !IGRAPH_EIT_END(eit);
    }

    /**
     * @brief Advance the iterator to the next edge.
     */
    void next() {
        IGRAPH_EIT_NEXT(eit);
    }

    /**
     * @brief Get the current edge.
     */
    void get(igraph_integer_t& from, igraph_integer_t& to) const {
        if (!has_next()) {
            throw std::out_of_range("EdgeIterator: No more edges to get");
        }
        igraph_integer_t edge_id = IGRAPH_EIT_GET(eit);
        igraph_edge(&igraph, edge_id, &from, &to);

        from = graph.get_original_node_id(from);
        to = graph.get_original_node_id(to);
    }
};
