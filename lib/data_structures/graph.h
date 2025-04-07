#pragma once

#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>

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

    // Compute the induced subgraph on a set of nodes
    Graph induced_subgraph(const std::vector<T>& nodes) const {
        Graph subgraph;
        for (const T& node : nodes) {
            subgraph.add_node(node);
        }
        for (const T& node : nodes) {
            for (const T& neighbor : get_neighbors(node)) {
                if (subgraph.has_node(neighbor)) {
                    subgraph.add_edge(node, neighbor);
                }
            }
        }
        return subgraph;
    }

    // Remove parallel edges and self-loops
    void remove_parallel_edges() {
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            std::unordered_map<int, bool> seen;
            auto& edges = adjacency_list[i];
            edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const Edge& edge) {
                if (seen[edge.to]) {
                    return true; // Remove this edge
                }
                seen[edge.to] = true;
                return false; // Keep this edge
            }), edges.end());
        }
    }

    // Remove self-loops
    void remove_self_loops() {
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            auto& edges = adjacency_list[i];
            edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const Edge& edge) {
                return edge.to == static_cast<int>(i); // Remove self-loop
            }), edges.end());
        }
    }

    // Remove floating nodes (nodes with no edges)
    void remove_floating_nodes() {
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            if (adjacency_list[i].empty()) {
                // Remove node
                int node_id = nodes[i];
                node_to_index.erase(node_id);
                nodes.erase(nodes.begin() + i);
                adjacency_list.erase(adjacency_list.begin() + i);
                --i; // Adjust index after removal
            }
        }
    }

    // Node iterator implementation
    class NodeIterator {
    private:
        typename std::vector<T>::const_iterator it;

    public:
        NodeIterator(typename std::vector<T>::const_iterator it) : it(it) {}
        
        const T& operator*() const { return *it; }
        NodeIterator& operator++() { ++it; return *this; }
        NodeIterator operator++(int) { NodeIterator tmp = *this; ++it; return tmp; }
        bool operator==(const NodeIterator& other) const { return it == other.it; }
        bool operator!=(const NodeIterator& other) const { return it != other.it; }
    };

    // Methods to get iterators
    NodeIterator begin() const { return NodeIterator(nodes.begin()); }
    NodeIterator end() const { return NodeIterator(nodes.end()); }
};
