#pragma once

#include <vector>
#include <algorithm> // For find, remove
#include <memory>
#include "graph.h" // Include your Graph header

/**
 * @brief A class that represents a subgraph pointer.
 * 
 * This class allows for the creation and manipulation of a subgraph
 * that is a subset of a larger graph. It provides methods to add and
 * remove vertices, check for edges, and retrieve neighbors.
 */
class SubgraphPointer {
private:
    Graph* main_graph_; // Pointer to the main graph
    std::unordered_set<int> vertices_; // Vertices in the subgraph, using int type

public:
    /**
     * @brief Constructor for SubgraphPointer.
     * 
     * @param main_graph Pointer to the main graph.
     * @param vertices Set of vertices in the subgraph.
     */
    SubgraphPointer(Graph* main_graph, const std::unordered_set<int>& vertices)
        : main_graph_(main_graph), vertices_(vertices) {}

    /**
     * @brief Check if a vertex is in the subgraph.
     * 
     * @param v The vertex ID to check.
     * @return True if the vertex is in the subgraph, false otherwise.
     */
    bool has_vertex(int v) const {
        return vertices_.find(v) != vertices_.end();
    }

    /**
     * @brief Check if the subgraph is empty.
     * 
     * @return True if the subgraph is empty, false otherwise.
     */
    size_t num_vertices() const {
        return vertices_.size();
    }

    /**
     * @brief Get the number of edges in the subgraph.
     * 
     * @return The number of edges in the subgraph.
     */
    size_t num_edges() const {
        size_t count = 0;
        for (const auto& v : vertices_) {
            count += get_neighbors(v).size();
        }
        return count / 2; // Each edge is counted twice
    }

    /** 
     * @brief Get the vertices in the subgraph.
     * 
     * @return A vector of vertex IDs in the subgraph.
     */
    std::unordered_set<int> get_vertices() const {
        return vertices_;
    }

    /** 
     * @brief Add a vertex to the subgraph.
     *
     * @param v The vertex ID to add.
     */
    void add_vertex(int v) {
        if (main_graph_->has_vertex(v) && !has_vertex(v)) {
            vertices_.push_back(v);
        }
    }

    /**
     * @brief Remove a vertex from the subgraph.
     *
     * @param v The vertex ID to remove.
     */
    void remove_vertex(int v) {
        vertices_.erase(std::remove(vertices_.begin(), vertices_.end(), v), vertices_.end());
    }

    /**
     * @brief Get the neighbors of a vertex in the subgraph.
     *
     * @param v The vertex ID whose neighbors to retrieve.
     * @return A vector of neighbor vertex IDs.
     */
    std::vector<int> get_neighbors(int v) const {
        if (!has_vertex(v)) return {};
        
        std::vector<int> neighbors;
        for (const auto& neighbor : main_graph_->get_neighbors(v)) {
            if (has_vertex(neighbor)) {
                neighbors.push_back(neighbor);
            }
        }
        return neighbors;
    }

    /**
     * @brief Check if an edge exists between two vertices in the subgraph.
     *
     * @param src The source vertex ID.
     * @param dst The destination vertex ID.
     * @return True if the edge exists, false otherwise.
     */
    bool has_edge(int src, int dst) const {
        return has_vertex(src) && has_vertex(dst) && main_graph_->has_edge(src, dst);
    }

    /**
     * @brief Add an edge between two vertices in the subgraph.
     *
     * @param src The source vertex ID.
     * @param dst The destination vertex ID.
     */
    void add_edge(int src, int dst) {
        add_vertex(src);
        add_vertex(dst);
        main_graph_->add_edge(src, dst);
    }

    // Remove edge operation is not needed for reference implementation
    // since edges are maintained by the main graph

    /**
     * @brief Get the degree of a vertex in the subgraph.
     *
     * @param v The vertex ID whose degree to retrieve.
     * @return The degree of the vertex.
     */
    size_t get_degree(int v) const {
        if (!has_vertex(v)) return 0;
        return get_neighbors(v).size();
    }

    /**
     * @brief Get the main graph.
     *
     * @return Pointer to the main graph.
     */
    Graph* get_main_graph() const {
        return main_graph_;
    }

    /**
     * @brief Get the subgraph as a new Graph object.
     *
     * @return A new Graph object representing the subgraph.
     */
    Graph realize() const {
        Graph subgraph;
        for (const auto& v : vertices_) {
            subgraph.add_vertex(v);
            for (const auto& neighbor : get_neighbors(v)) {
                subgraph.add_edge(v, neighbor);
            }
        }
        return subgraph;
    }
};
