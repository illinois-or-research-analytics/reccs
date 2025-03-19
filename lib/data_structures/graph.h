#pragma once

#include <vector>
#include <unordered_map>
#include <stdexcept>

template <typename T>
class Graph {
private:
    struct Edge {
        int to;
        
        Edge(int to) : to(to) {}
    };
    
    std::unordered_map<T, int> node_to_index;
    std::vector<T> nodes;
    std::vector<std::vector<Edge>> adjacency_list;

public:
    Graph() = default;
    
    // Add a node to the graph
    void add_node(const T& node) {
        if (node_to_index.find(node) == node_to_index.end()) {
            node_to_index[node] = nodes.size();
            nodes.push_back(node);
            adjacency_list.push_back(std::vector<Edge>());
        }
    }
    
    // Add an edge between two nodes
    void add_edge(const T& from, const T& to) {
        // Add nodes if they don't exist
        add_node(from);
        add_node(to);
        
        int from_idx = node_to_index[from];
        int to_idx = node_to_index[to];
        
        // Add edge
        adjacency_list[from_idx].push_back(Edge(to_idx));
    }
    
    // Check if a node exists
    bool has_node(const T& node) const {
        return node_to_index.find(node) != node_to_index.end();
    }
    
    // Get all nodes
    const std::vector<T>& get_nodes() const {
        return nodes;
    }
    
    // Get all neighbors of a node
    std::vector<T> get_neighbors(const T& node) const {
        if (!has_node(node)) {
            throw std::runtime_error("Node not found in graph");
        }
        
        int node_idx = node_to_index.at(node);
        std::vector<T> neighbors;
        
        for (const Edge& edge : adjacency_list[node_idx]) {
            neighbors.push_back(nodes[edge.to]);
        }
        
        return neighbors;
    }

    // Get the number of edges in the graph
    size_t edge_count() const {
        size_t count = 0;
        for (const auto& edges : adjacency_list) {
            count += edges.size();
        }
        return count;
    }

    // Get the number of nodes in the graph
    size_t node_count() const {
        return nodes.size();
    }
};
