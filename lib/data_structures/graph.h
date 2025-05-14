#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <cassert>
#include <iostream>

/**
 * @class Graph
 * @brief Implementation of the Double-Index (DI) graph data structure for undirected, unweighted graphs
 * 
 * This class implements the Double-Index data structure described in the paper
 * "Interactive Graph Stream Analytics in Arkouda". It provides efficient lookup
 * from edges to vertices and from vertices to their adjacency lists.
 */
class Graph {
public:
    /**
     * @brief Default constructor
     */
    Graph() : n_vertices_(0), n_edges_(0) {}
    
    /**
     * @brief Constructor with parameters
     * 
     * @param n_vertices Number of vertices
     * @param n_edges Number of edges
     */
    Graph(size_t n_vertices, size_t n_edges) :
        n_vertices_(n_vertices),
        n_edges_(n_edges) {
        
        // Initialize arrays with appropriate sizes
        src_.resize(n_edges);
        dst_.resize(n_edges);
        
        start_i_.resize(n_vertices, -1);  // -1 indicates no outgoing edges
        neighbor_.resize(n_vertices, 0);
        
        // For undirected graphs, initialize reverse arrays
        src_r_.resize(n_edges);
        dst_r_.resize(n_edges);
        start_i_r_.resize(n_vertices, -1);
        neighbor_r_.resize(n_vertices, 0);
    }
    
    /**
     * @brief Builds the graph from edge lists with explicit node mapping
     * 
     * @param src Source vertex IDs (original non-continuous IDs)
     * @param dst Destination vertex IDs (original non-continuous IDs)
     * @param node_map Optional node mapping from original IDs to continuous 0-based indices
     */
    void build_from_edges(const std::vector<int>& src, 
                          const std::vector<int>& dst,
                          const std::unordered_map<int, int>& node_map = std::unordered_map<int, int>()) {
        assert(src.size() == dst.size() && "Source and destination arrays must have the same size");
        
        n_edges_ = src.size();
        
        // If node mapping is provided, use it, otherwise build one
        if (!node_map.empty()) {
            node_mapping_ = node_map;
            reverse_mapping_.resize(node_map.size());
            
            for (const auto& pair : node_map) {
                reverse_mapping_[pair.second] = pair.first;
            }
            
            n_vertices_ = node_map.size();
        } else {
            // Collect all unique vertex IDs
            std::unordered_set<int> vertices;
            for (size_t i = 0; i < src.size(); i++) {
                vertices.insert(src[i]);
                vertices.insert(dst[i]);
            }
            
            // Create mapping from original ID to continuous 0-based index
            int idx = 0;
            for (int vertex : vertices) {
                node_mapping_[vertex] = idx;
                reverse_mapping_.push_back(vertex);
                idx++;
            }
            
            n_vertices_ = vertices.size();
        }
        
        // Map original IDs to continuous indices
        std::vector<int> mapped_src(n_edges_);
        std::vector<int> mapped_dst(n_edges_);
        
        #pragma omp parallel for
        for (size_t i = 0; i < n_edges_; i++) {
            mapped_src[i] = node_mapping_[src[i]];
            mapped_dst[i] = node_mapping_[dst[i]];
        }
        
        src_ = mapped_src;
        dst_ = mapped_dst;
        
        // Initialize vertex arrays
        start_i_.resize(n_vertices_, -1);
        neighbor_.resize(n_vertices_, 0);
        
        // Sort edges by source vertex ID for building the vertex index
        std::vector<size_t> indices(n_edges_);
        #pragma omp parallel for
        for (size_t i = 0; i < n_edges_; i++) {
            indices[i] = i;
        }
        
        std::sort(indices.begin(), indices.end(), [this](size_t i, size_t j) {
            return (src_[i] < src_[j]) || (src_[i] == src_[j] && dst_[i] < dst_[j]);
        });
        
        // Reorder the edge arrays according to sorted indices
        std::vector<int> sorted_src(n_edges_);
        std::vector<int> sorted_dst(n_edges_);
        
        #pragma omp parallel for
        for (size_t i = 0; i < n_edges_; i++) {
            sorted_src[i] = src_[indices[i]];
            sorted_dst[i] = dst_[indices[i]];
        }
        
        src_ = sorted_src;
        dst_ = sorted_dst;
        
        // Build the vertex index (STR and NEI arrays from the paper)
        for (size_t i = 0; i < n_edges_; i++) {
            int source = src_[i];
            if (start_i_[source] == -1) {
                start_i_[source] = i;
            }
            neighbor_[source]++;
        }
        
        // For undirected graphs, build the reverse arrays
        build_reverse_arrays();
    }
    
    /**
     * @brief Build the reverse arrays for undirected access
     */
    void build_reverse_arrays() {
        src_r_ = dst_;  // Reverse the source and destination
        dst_r_ = src_;
        
        start_i_r_.resize(n_vertices_, -1);
        neighbor_r_.resize(n_vertices_, 0);
        
        // Sort by the new source (original destination)
        std::vector<size_t> indices(n_edges_);
        #pragma omp parallel for
        for (size_t i = 0; i < n_edges_; i++) {
            indices[i] = i;
        }
        
        std::sort(indices.begin(), indices.end(), [this](size_t i, size_t j) {
            return (src_r_[i] < src_r_[j]) || (src_r_[i] == src_r_[j] && dst_r_[i] < dst_r_[j]);
        });
        
        // Reorder the reverse edge arrays
        std::vector<int> sorted_src_r(n_edges_);
        std::vector<int> sorted_dst_r(n_edges_);
        
        #pragma omp parallel for
        for (size_t i = 0; i < n_edges_; i++) {
            sorted_src_r[i] = src_r_[indices[i]];
            sorted_dst_r[i] = dst_r_[indices[i]];
        }
        
        src_r_ = sorted_src_r;
        dst_r_ = sorted_dst_r;
        
        // Build the reverse vertex index
        for (size_t i = 0; i < n_edges_; i++) {
            int source = src_r_[i];
            if (start_i_r_[source] == -1) {
                start_i_r_[source] = i;
            }
            neighbor_r_[source]++;
        }
    }
    
    /**
     * @brief Gets the neighbors of a vertex (using internal continuous ID)
     * 
     * @param vertex_id The internal vertex ID (0-based continuous index)
     * @return std::vector<int> List of neighbor vertex IDs (internal continuous IDs)
     */
    std::vector<int> get_neighbors(int vertex_id) const {
        std::vector<int> neighbors;
        
        // Get outgoing edges
        if (start_i_[vertex_id] != -1) {
            size_t start = start_i_[vertex_id];
            size_t count = neighbor_[vertex_id];
            
            neighbors.reserve(count);
            for (size_t i = 0; i < count; i++) {
                neighbors.push_back(dst_[start + i]);
            }
        }
        
        // Get incoming edges
        if (start_i_r_[vertex_id] != -1) {
            size_t start = start_i_r_[vertex_id];
            size_t count = neighbor_r_[vertex_id];
            
            neighbors.reserve(neighbors.size() + count);
            for (size_t i = 0; i < count; i++) {
                neighbors.push_back(dst_r_[start + i]);
            }
        }
        
        return neighbors;
    }
    
    /**
     * @brief Gets the neighbors of a vertex using original ID
     * 
     * @param original_id The original vertex ID
     * @return std::vector<int> List of neighbor vertex IDs (original IDs)
     */
    std::vector<int> get_neighbors_original(int original_id) const {
        if (node_mapping_.find(original_id) == node_mapping_.end()) {
            return {};
        }
        
        int internal_id = node_mapping_.at(original_id);
        std::vector<int> internal_neighbors = get_neighbors(internal_id);
        
        // Convert internal IDs back to original IDs
        std::vector<int> original_neighbors;
        original_neighbors.reserve(internal_neighbors.size());
        
        for (int neighbor : internal_neighbors) {
            original_neighbors.push_back(reverse_mapping_[neighbor]);
        }
        
        return original_neighbors;
    }
    
    /**
     * @brief Gets the edge ID between two vertices if it exists
     * 
     * @param source Source vertex ID (internal ID)
     * @param target Target vertex ID (internal ID)
     * @return int Edge ID or -1 if no edge exists
     */
    int get_edge_id(int source, int target) const {
        if (start_i_[source] == -1) return -1;
        
        size_t start = start_i_[source];
        size_t count = neighbor_[source];
        
        for (size_t i = 0; i < count; i++) {
            if (dst_[start + i] == target) {
                return start + i;
            }
        }
        
        return -1;
    }
    
    /**
     * @brief Gets the edge ID between two vertices using original IDs
     * 
     * @param original_source Original source vertex ID
     * @param original_target Original target vertex ID
     * @return int Edge ID or -1 if no edge exists
     */
    int get_edge_id_original(int original_source, int original_target) const {
        if (node_mapping_.find(original_source) == node_mapping_.end() ||
            node_mapping_.find(original_target) == node_mapping_.end()) {
            return -1;
        }
        
        int internal_source = node_mapping_.at(original_source);
        int internal_target = node_mapping_.at(original_target);
        
        return get_edge_id(internal_source, internal_target);
    }
    
    /**
     * @brief Convert internal vertex ID to original ID
     * 
     * @param internal_id Internal vertex ID (0-based continuous index)
     * @return int Original vertex ID
     */
    int to_original_id(int internal_id) const {
        if (internal_id < 0 || internal_id >= static_cast<int>(n_vertices_)) {
            return -1;  // Invalid ID
        }
        return reverse_mapping_[internal_id];
    }
    
    /**
     * @brief Convert original vertex ID to internal ID
     * 
     * @param original_id Original vertex ID
     * @return int Internal vertex ID (0-based continuous index), or -1 if not found
     */
    int to_internal_id(int original_id) const {
        auto it = node_mapping_.find(original_id);
        if (it == node_mapping_.end()) {
            return -1;  // Not found
        }
        return it->second;
    }
    
    /**
     * @brief Print graph statistics
     */
    void print_stats() const {
        std::cout << "Graph Statistics:" << std::endl;
        std::cout << "  Vertices: " << n_vertices_ << std::endl;
        std::cout << "  Edges: " << n_edges_ << std::endl;
        std::cout << "  Has node mapping: " << (!node_mapping_.empty() ? "Yes" : "No") << std::endl;
    }
    
    // Accessors
    size_t num_vertices() const { return n_vertices_; }
    size_t num_edges() const { return n_edges_; }
    
    const std::vector<int>& src() const { return src_; }
    const std::vector<int>& dst() const { return dst_; }
    const std::vector<int>& start_i() const { return start_i_; }
    const std::vector<int>& neighbor() const { return neighbor_; }
    
    const std::vector<int>& src_r() const { return src_r_; }
    const std::vector<int>& dst_r() const { return dst_r_; }
    const std::vector<int>& start_i_r() const { return start_i_r_; }
    const std::vector<int>& neighbor_r() const { return neighbor_r_; }
    
    const std::unordered_map<int, int>& node_mapping() const { return node_mapping_; }
    const std::vector<int>& reverse_mapping() const { return reverse_mapping_; }
    
private:
    size_t n_vertices_;    // Number of vertices
    size_t n_edges_;       // Number of edges
    
    // Edge index arrays
    std::vector<int> src_;          // Source vertex of each edge
    std::vector<int> dst_;          // Destination vertex of each edge
    
    // Vertex index arrays
    std::vector<int> start_i_;      // Starting index of edges for each vertex
    std::vector<int> neighbor_;     // Number of neighbors for each vertex
    
    // Reverse arrays for undirected access
    std::vector<int> src_r_;        // Reverse source vertex of each edge
    std::vector<int> dst_r_;        // Reverse destination vertex of each edge
    std::vector<int> start_i_r_;    // Reverse starting index
    std::vector<int> neighbor_r_;   // Reverse number of neighbors
    
    // Node mapping
    std::unordered_map<int, int> node_mapping_;  // Original ID -> Internal ID
    std::vector<int> reverse_mapping_;           // Internal ID -> Original ID
};
