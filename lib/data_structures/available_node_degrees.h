#ifndef AVAILABLE_NODE_DEGREES_H
#define AVAILABLE_NODE_DEGREES_H

#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <atomic>
#include <memory>
#include "graph.h"

/**
 * Thread-safe manager for tracking available node degrees across all clusters
 * This replicates the Python logic where nodes have "budget" for additional edges
 */
class AvailableNodeDegreesManager {
private:
    // Maps node_id -> available degree budget
    std::unordered_map<uint64_t, std::atomic<int32_t>> available_degrees;
    
    // Fast lookup set of nodes with available degrees > 0
    std::unordered_set<uint64_t> available_nodes;
    mutable std::mutex available_nodes_mutex;
    
    // Reference graphs for degree computation
    std::shared_ptr<Graph> empirical_graph;
    std::shared_ptr<Graph> synthetic_graph;
    
    // Node mapping between graphs
    std::unordered_map<uint64_t, uint32_t> emp_node_mapping;
    std::unordered_map<uint64_t, uint32_t> syn_node_mapping;

public:
    AvailableNodeDegreesManager(
        std::shared_ptr<Graph> emp_graph,
        std::shared_ptr<Graph> syn_graph,
        const std::unordered_map<uint64_t, uint32_t>& emp_mapping,
        const std::unordered_map<uint64_t, uint32_t>& syn_mapping)
        : empirical_graph(emp_graph), synthetic_graph(syn_graph),
          emp_node_mapping(emp_mapping), syn_node_mapping(syn_mapping) {
        initialize_available_degrees();
    }
    
    /**
     * Initialize available degrees based on empirical vs synthetic degree differences
     * Replicates: deg_diff = emp_degrees[c_node] - sbm_degrees[c_node]
     */
    void initialize_available_degrees() {
        std::lock_guard<std::mutex> lock(available_nodes_mutex);
        
        for (const auto& [node_id, emp_node] : emp_node_mapping) {
            auto syn_it = syn_node_mapping.find(node_id);
            if (syn_it == syn_node_mapping.end()) continue;
            
            uint32_t syn_node = syn_it->second;
            
            // Get degrees from both graphs
            uint32_t emp_degree = get_degree(*empirical_graph, emp_node);
            uint32_t syn_degree = get_degree(*synthetic_graph, syn_node);
            
            // Calculate available degree budget
            int32_t deg_diff = static_cast<int32_t>(emp_degree) - static_cast<int32_t>(syn_degree);
            
            if (deg_diff > 0) {
                available_degrees[node_id].store(deg_diff);
                available_nodes.insert(node_id);
            }
        }
        
        std::cout << "Initialized " << available_nodes.size() 
                  << " nodes with available degrees" << std::endl;
    }
    
    /**
     * Try to consume available degree budget for a node
     * Returns true if successful, false if insufficient budget
     */
    bool try_consume_degree(uint64_t node_id, int32_t amount = 1) {
        auto it = available_degrees.find(node_id);
        if (it == available_degrees.end()) return false;
        
        // Atomic compare-and-swap to consume budget
        int32_t current = it->second.load();
        while (current >= amount) {
            if (it->second.compare_exchange_weak(current, current - amount)) {
                // Successfully consumed, check if node should be removed from available set
                if (current - amount <= 0) {
                    std::lock_guard<std::mutex> lock(available_nodes_mutex);
                    available_nodes.erase(node_id);
                }
                return true;
            }
        }
        return false;
    }
    
    /**
     * Get available degree budget for a node
     */
    int32_t get_available_degree(uint64_t node_id) const {
        auto it = available_degrees.find(node_id);
        return (it != available_degrees.end()) ? it->second.load() : 0;
    }
    
    /**
     * Get set of nodes with available degrees (thread-safe copy)
     */
    std::unordered_set<uint64_t> get_available_nodes() const {
        std::lock_guard<std::mutex> lock(available_nodes_mutex);
        return available_nodes;  // Copy
    }
    
    /**
     * Get intersection of available nodes with a given node set
     * Replicates: available_c_nodes = available_node_set.intersection(nodes)
     */
    std::unordered_set<uint64_t> get_available_intersection(
        const std::unordered_set<uint64_t>& cluster_nodes) const {
        
        std::lock_guard<std::mutex> lock(available_nodes_mutex);
        std::unordered_set<uint64_t> result;
        
        for (uint64_t node : cluster_nodes) {
            if (available_nodes.count(node) > 0) {
                result.insert(node);
            }
        }
        return result;
    }
    
    /**
     * Remove nodes from available set without consuming degrees
     * (for when nodes are exhausted through other means)
     */
    void remove_from_available(const std::unordered_set<uint64_t>& nodes_to_remove) {
        std::lock_guard<std::mutex> lock(available_nodes_mutex);
        for (uint64_t node : nodes_to_remove) {
            available_nodes.erase(node);
        }
    }
    
    /**
     * Get statistics about available degrees
     */
    struct AvailableStats {
        size_t total_available_nodes;
        int32_t total_available_degrees;
        double avg_available_degree;
    };
    
    AvailableStats get_stats() const {
        std::lock_guard<std::mutex> lock(available_nodes_mutex);
        AvailableStats stats;
        stats.total_available_nodes = available_nodes.size();
        stats.total_available_degrees = 0;
        
        for (uint64_t node : available_nodes) {
            auto it = available_degrees.find(node);
            if (it != available_degrees.end()) {
                stats.total_available_degrees += it->second.load();
            }
        }
        
        stats.avg_available_degree = (stats.total_available_nodes > 0) ? 
            static_cast<double>(stats.total_available_degrees) / stats.total_available_nodes : 0.0;
            
        return stats;
    }

private:
    uint32_t get_degree(const Graph& g, uint32_t node) const {
        return g.row_ptr[node + 1] - g.row_ptr[node];
    }
};

/**
 * Modified GraphTask to include reference to available degrees manager
 */
struct GraphTaskWithDegrees {
    // Original task data
    std::shared_ptr<Graph> subgraph;
    TaskType task_type;
    std::string cluster_id;
    uint32_t cluster_idx;
    uint32_t min_degree_requirement;
    
    // New: reference to shared degree manager
    std::shared_ptr<AvailableNodeDegreesManager> degree_manager;
    
    // Cluster node IDs (for intersection operations)
    std::unordered_set<uint64_t> cluster_node_ids;
    
    GraphTaskWithDegrees(
        std::shared_ptr<Graph> g, 
        TaskType t, 
        const std::string& cid,
        uint32_t cidx, 
        uint32_t min_deg,
        std::shared_ptr<AvailableNodeDegreesManager> deg_mgr,
        const std::unordered_set<uint64_t>& node_ids)
        : subgraph(g), task_type(t), cluster_id(cid), cluster_idx(cidx),
          min_degree_requirement(min_deg), degree_manager(deg_mgr),
          cluster_node_ids(node_ids) {}
};

/**
 * Helper functions for stage implementations
 */
namespace degree_aware_stages {

/**
 * Get available nodes for this cluster (intersection operation)
 */
std::unordered_set<uint64_t> get_cluster_available_nodes(
    const GraphTaskWithDegrees& task) {
    return task.degree_manager->get_available_intersection(task.cluster_node_ids);
}

/**
 * Try to add an edge with degree budget checking
 * Returns true if edge was added and budgets consumed
 */
bool try_add_edge_with_budget(
    Graph& g, 
    uint32_t u, uint32_t v,
    const GraphTaskWithDegrees& task) {
    
    // Convert subgraph indices back to original node IDs
    uint64_t u_id = g.id_map[u];
    uint64_t v_id = g.id_map[v];
    
    // Check if both nodes have available degree budget
    bool u_available = task.degree_manager->get_available_degree(u_id) > 0;
    bool v_available = task.degree_manager->get_available_degree(v_id) > 0;
    
    // Add the edge
    std::vector<std::pair<uint32_t, uint32_t>> edges = {{u, v}};
    add_edges_batch(g, edges);
    
    // Consume budgets if available
    bool degree_corrected = false;
    if (u_available && task.degree_manager->try_consume_degree(u_id)) {
        degree_corrected = true;
    }
    if (v_available && task.degree_manager->try_consume_degree(v_id)) {
        degree_corrected = true;
    }
    
    return degree_corrected;
}

/**
 * Select edge endpoint with preference for available degree nodes
 */
uint32_t select_edge_endpoint(
    const std::vector<uint32_t>& candidates,
    const Graph& g,
    const GraphTaskWithDegrees& task,
    std::mt19937& gen) {
    
    if (candidates.empty()) return UINT32_MAX;
    
    // Separate candidates into those with/without available degrees
    std::vector<uint32_t> available_candidates;
    std::vector<uint32_t> other_candidates;
    
    for (uint32_t candidate : candidates) {
        uint64_t candidate_id = g.id_map[candidate];
        if (task.degree_manager->get_available_degree(candidate_id) > 0) {
            available_candidates.push_back(candidate);
        } else {
            other_candidates.push_back(candidate);
        }
    }
    
    // Prefer available degree candidates
    if (!available_candidates.empty()) {
        std::uniform_int_distribution<size_t> dist(0, available_candidates.size() - 1);
        return available_candidates[dist(gen)];
    } else if (!other_candidates.empty()) {
        std::uniform_int_distribution<size_t> dist(0, other_candidates.size() - 1);
        return other_candidates[dist(gen)];
    }
    
    return UINT32_MAX;
}

} // namespace degree_aware_stages

#endif // AVAILABLE_NODE_DEGREES_H
