#pragma once

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>

#include <igraph/igraph.h>


/**
 * @brief Wrapper class for igraph_t to represent a graph.
 */
class Graph {
private:
    igraph_t graph;  // The underlying igraph object
    int num_nodes;  // Number of nodes in the graph
    int num_edges;  // Number of edges in the graph

public:
    /**
     * @brief Default constructor.
     */
    Graph() : num_nodes(0), num_edges(0) {
        igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    }

    /**
     * @brief Constructor that initializes the graph from an igraph_t object.
     */
    Graph(const igraph_t& g) {
        // Initialize our graph with a copy
        igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
        igraph_copy(&graph, &g);
        
        // Now get the counts
        num_nodes = igraph_vcount(&graph);
        num_edges = igraph_ecount(&graph);
        
        std::cout << "Graph constructor called" << std::endl;
        std::cout << "Number of nodes: " << num_nodes << std::endl;
        std::cout << "Number of edges: " << num_edges << std::endl;
    }

    /**
     * @brief Get the number of nodes in the graph.
     * 
     * @return The number of nodes.
     */
    int get_num_nodes() const {
        return num_nodes;
    }

    /**
     * @brief Get the number of edges in the graph.
     * 
     * @return The number of edges.
     */
    int get_num_edges() const {
        return num_edges;
    }   

    /**
     * @brief Get the underlying igraph object.
     * 
     * @return A reference to the igraph object.
     */
    const igraph_t& get_igraph() const {
        return graph;
    }

    /**
     * @brief Get the degree of a node.
     * 
     * @param node_id The ID of the node.
     * @return The degree of the node.
     */
    int get_degree(int node_id) const {
        igraph_vector_int_t degree;
        igraph_vector_int_init(&degree, 1);
        igraph_degree(&graph, &degree, igraph_vss_1(node_id), IGRAPH_ALL, 0);
        int result = VECTOR(degree)[0];
        igraph_vector_int_destroy(&degree);
        return result;
    }

    /**
     * @brief Remove parallel edges and self-loops from the graph.
     * 
     * @param remove_self_loops Whether to remove self-loops.
     * @param remove_parallel_edges Whether to remove parallel edges.
     */
    void cleanup(bool remove_self_loops = true, bool remove_parallel_edges = true) {
        igraph_bool_t multiple = remove_parallel_edges;
        igraph_bool_t loops = remove_self_loops;
        igraph_attribute_combination_t comb;
        igraph_attribute_combination_init(&comb);
        igraph_simplify(&graph, multiple, loops, &comb);
        igraph_attribute_combination_destroy(&comb);
        
        // Update the counts after simplification
        num_nodes = igraph_vcount(&graph);
        num_edges = igraph_ecount(&graph);
    }

    /**
     * @brief Get the induced subgraph of a set of vertices.
     *
     * @param vertices A vector of vertex IDs.
     * @return A new Graph object representing the induced subgraph.
     */
    Graph get_induced_subgraph(const std::vector<int>& vertices) const {
        igraph_t subgraph;
        igraph_vector_int_t vertex_vector;
        igraph_vs_t vs;
        
        igraph_vector_int_init(&vertex_vector, vertices.size());
        
        for (size_t i = 0; i < vertices.size(); ++i) {
            VECTOR(vertex_vector)[i] = vertices[i];
        }
        
        igraph_vs_vector(&vs, &vertex_vector);
        igraph_induced_subgraph(&graph, &subgraph, vs, IGRAPH_SUBGRAPH_AUTO);
        
        Graph result(subgraph);
        igraph_vs_destroy(&vs);
        igraph_vector_int_destroy(&vertex_vector);
        
        return result;
    }

    /**
     * @brief Get the original node ID for a node in the graph.
     * 
     * @param node_id The current ID of the node.
     * @return The original ID of the node.
     */
    int get_original_node_id(int node_id) const {
        // For now, just return the same ID
        // This will be implemented properly later
        return node_id;
    }

    /**
     * @brief Set the mapping between current and original node IDs.
     * 
     * @param mapping The mapping from current to original node IDs.
     */
    void set_node_id_mapping(const std::unordered_map<int, int>& mapping) {
        // This will be implemented later
    }
};
