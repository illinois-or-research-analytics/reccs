#pragma once

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <limits>

#include <igraph/igraph.h>


/**
 * @brief Wrapper class for igraph_t to represent a graph.
 */
class Graph {
private:
    igraph_t graph;  // The underlying igraph object
    std::unordered_map<int, int> id_to_index;  // Maps original node IDs to contiguous indices
    std::unordered_map<int, int> index_to_id;  // Maps contiguous indices back to original node IDs
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
    * @brief Contstructor from an igraph object.
     */
    Graph(const igraph_t& g, std::unordered_map<int, int> id_to_index) : graph(g), id_to_index(id_to_index) {
        num_nodes = igraph_vcount(&graph);
        num_edges = igraph_ecount(&graph);

        // Use id_to_index to compute index_to_id
        index_to_id.reserve(num_nodes);
        for (const auto& pair : id_to_index) {
            index_to_id[pair.second] = pair.first;
        }
    }

    /**
     * @brief Destructor.
     */
    ~Graph() {
        igraph_destroy(&graph);
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
     * @brief Get the original node ID for a given index.
     * 
     * @param index The index of the node.
     * @return The original node ID.
     * @throws std::out_of_range if the index is out of bounds.
     */
    int get_original_node_id(int index) const {
        if (index < 0 || index >= num_nodes) {
            throw std::out_of_range("Index out of bounds");
        }
        return index_to_id.at(index);
    }

    /**
     * @brief Get the index for a given original node ID.
     * 
     * @param id The original node ID.
     * @return The index of the node.
     * @throws std::out_of_range if the ID is not found.
     */
    int get_index(int id) const {
        auto it = id_to_index.find(id);
        if (it == id_to_index.end()) {
            throw std::out_of_range("ID not found");
        }
        return it->second;
    }

    /**
     * @brief Add an edge to the graph.
     * 
     * @param from The ID of the source node.
     * @param to The ID of the target node.
     */
    void add_edge(int from, int to) {
        if (id_to_index.find(from) == id_to_index.end()) {
            id_to_index[from] = num_nodes++;
            index_to_id[num_nodes - 1] = from;
        }
        if (id_to_index.find(to) == id_to_index.end()) {
            id_to_index[to] = num_nodes++;
            index_to_id[num_nodes - 1] = to;
        }
        igraph_add_edge(&graph, id_to_index[from], id_to_index[to]);
        num_edges++;
    }

    /**
     * @brief Get the underlying igraph object.
     * 
     * @return A reference to the igraph object.
     */
    igraph_t& get_igraph() {
        return graph;
    }

    /**
     * @brief Get the underlying igraph object (const version).
     * 
     * @return A const reference to the igraph object.
     */
    const igraph_t& get_igraph() const {
        return graph;
    }

    /**
     * @brief Get the degree of a node.
     * 
     * @param id The ID of the node.
     * @return The degree of the node.
     */
    int get_degree(int id) const {
        int index = get_index(id);
        igraph_vector_int_t degrees;
        igraph_vector_int_init(&degrees, 1);
        igraph_degree(&graph, &degrees, igraph_vss_1(index), IGRAPH_ALL, IGRAPH_NO_LOOPS);
        int result = VECTOR(degrees)[0];
        igraph_vector_int_destroy(&degrees);
        return result;
    }

    /**
     * @brief Get the neighbors of a node.
     * 
     * @param id The ID of the node.
     * @return A vector of neighbor IDs.
     */
    std::vector<int> get_neighbors(int id) const {
        int index = get_index(id);
        igraph_vector_int_t neighbors;
        igraph_vector_int_init(&neighbors, 0);
        igraph_neighbors(&graph, &neighbors, index, IGRAPH_ALL);
        
        std::vector<int> result;
        for (int i = 0; i < igraph_vector_int_size(&neighbors); ++i) {
            result.push_back(index_to_id.at(VECTOR(neighbors)[i]));
        }
        
        igraph_vector_int_destroy(&neighbors);
        return result;
    }

    /**
     * @brief Get the induced subgraph from a set of vertices.
     *
     * @param vertices A vector of vertex indices.
     * @return A new Graph object representing the induced subgraph.
    */
    Graph get_induced_subgraph(const std::vector<int>& vertices) const {
        igraph_vs_t vs;
        igraph_vector_int_t vertex_vector;
        igraph_vector_int_init(&vertex_vector, vertices.size());
        
        for (size_t i = 0; i < vertices.size(); ++i) {
            VECTOR(vertex_vector)[i] = get_index(vertices[i]);
        }
        
        igraph_vs_vector(&vs, &vertex_vector);
        
        igraph_t subgraph;
        igraph_induced_subgraph(&graph, &subgraph, vs, IGRAPH_SUBGRAPH_AUTO);
        
        Graph result;
        result.graph = subgraph;
        result.num_nodes = vertices.size();
        result.num_edges = igraph_ecount(&subgraph);
        
        // Map original IDs to new indices
        for (size_t i = 0; i < vertices.size(); ++i) {
            result.id_to_index[vertices[i]] = i;
            result.index_to_id[i] = vertices[i];
        }
        
        igraph_vs_destroy(&vs);
        igraph_vector_int_destroy(&vertex_vector);
        
        return result;
    }
};