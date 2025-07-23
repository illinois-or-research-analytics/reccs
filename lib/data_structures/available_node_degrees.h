#ifndef AVAILABLE_NODE_DEGREES_H
#define AVAILABLE_NODE_DEGREES_H

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fstream>
#include <iostream>
#include "graph.h"

using json = nlohmann::json;

/**
 * Global immutable degree mapping - loaded from JSON, never modified
 * Each task gets its own local copy of relevant degrees
 */
class AvailableNodeDegreesManager {
private:
    // Global immutable mapping: node_id -> available_degree_budget
    std::unordered_map<uint64_t, int32_t> global_degree_budgets;

public:
    /**
     * Constructor that loads degree deficits from JSON file
     */
    AvailableNodeDegreesManager(const std::string& json_filename) {
        load_degree_budgets_from_json(json_filename);
    }
    
    /**
     * Load degree budgets from JSON file - called once during initialization
     */
    void load_degree_budgets_from_json(const std::string& json_filename) {
        std::cout << "Loading degree budgets from JSON: " << json_filename << std::endl;
        
        std::ifstream file(json_filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open JSON file: " + json_filename);
        }
        
        json j;
        file >> j;
        file.close();
        
        size_t positive_budgets = 0;
        size_t nodes_loaded = 0;
        
        // Load node_id -> deficit mapping from JSON
        for (const auto& [node_id_str, deficit] : j.items()) {
            try {
                uint64_t node_id = std::stoull(node_id_str);
                int32_t budget = deficit.get<int32_t>();
                
                // Only store positive budgets (matching original Python logic)
                if (budget > 0) {
                    global_degree_budgets[node_id] = budget;
                    positive_budgets++;
                }
                
                nodes_loaded++;
                
                // Debug output for first few nodes
                if (nodes_loaded <= 10) {
                    std::cout << "  Node " << node_id << ": budget=" << budget << std::endl;
                }
                
            } catch (const std::exception& e) {
                std::cerr << "Warning: Failed to parse node " << node_id_str 
                          << " with deficit " << deficit << ": " << e.what() << std::endl;
            }
        }
        
        std::cout << "Loaded degree budgets: " << nodes_loaded << " total nodes processed, " 
                  << positive_budgets << " with positive budget stored" << std::endl;
    }
    
    /**
     * Create local degree budget map for a specific cluster
     * Each task gets its own independent copy - NO SHARED STATE!
     */
    std::unordered_map<uint64_t, int32_t> create_local_degrees(
        const std::unordered_set<uint64_t>& cluster_node_ids) const {
        
        std::unordered_map<uint64_t, int32_t> local_budgets;
        
        // Extract only the nodes belonging to this cluster that have positive budgets
        for (uint64_t node_id : cluster_node_ids) {
            auto it = global_degree_budgets.find(node_id);
            if (it != global_degree_budgets.end()) {
                local_budgets[node_id] = it->second;  // Copy the positive budget
            }
            // Note: nodes with zero or negative budgets are not stored
        }
        
        std::cout << "Created local degrees for cluster with " << local_budgets.size() 
                  << " nodes with positive budget" << std::endl;
        
        return local_budgets;  // Return by value - each task gets its own copy
    }
    
    /**
     * Get list of nodes with positive budget for a cluster
     */
    std::vector<uint64_t> get_available_nodes_for_cluster(
        const std::unordered_set<uint64_t>& cluster_node_ids) const {
        
        std::vector<uint64_t> available_nodes;
        
        for (uint64_t node_id : cluster_node_ids) {
            auto it = global_degree_budgets.find(node_id);
            if (it != global_degree_budgets.end()) {
                available_nodes.push_back(node_id);
            }
        }
        
        return available_nodes;
    }
    
    int32_t get_available_degree(uint64_t node_id) const {
        auto it = global_degree_budgets.find(node_id);
        return (it != global_degree_budgets.end()) ? it->second : 0;
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
        AvailableStats stats;
        stats.total_available_nodes = global_degree_budgets.size();
        stats.total_available_degrees = 0;
        
        for (const auto& [node_id, budget] : global_degree_budgets) {
            stats.total_available_degrees += budget;
        }
        
        stats.avg_available_degree = (stats.total_available_nodes > 0) ? 
            static_cast<double>(stats.total_available_degrees) / stats.total_available_nodes : 0.0;
            
        return stats;
    }
};

/**
 * Enhanced GraphTask with local degree management - NO SHARED STATE!
 */
struct GraphTaskWithDegrees {
    // Original task data
    std::shared_ptr<Graph> subgraph;
    TaskType task_type;
    std::string cluster_id;
    uint32_t cluster_idx;
    uint32_t min_degree_requirement;
    
    // Reference to global degree manager (read-only after initialization)
    std::shared_ptr<AvailableNodeDegreesManager> degree_manager;
    
    // Cluster node IDs
    std::unordered_set<uint64_t> cluster_node_ids;
    
    // LOCAL degree budgets - each task has its own copy!
    mutable std::unordered_map<uint64_t, int32_t> local_degree_budgets;
    mutable std::vector<uint64_t> local_available_nodes;
    mutable bool local_degrees_initialized = false;
    
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
    
    /**
     * Initialize local degrees for this task (called once per task)
     */
    void initialize_local_degrees() const {
        if (local_degrees_initialized) return;
        
        // Create local copy of degree budgets for this cluster
        local_degree_budgets = degree_manager->create_local_degrees(cluster_node_ids);
        
        // Create initial list of available nodes (only computed once!)
        local_available_nodes.reserve(local_degree_budgets.size());
        for (const auto& [node_id, budget] : local_degree_budgets) {
            if (budget > 0) {
                local_available_nodes.push_back(node_id);
            }
        }
        
        std::cout << "Initialized local degrees with " << local_available_nodes.size() 
                  << " initially available nodes" << std::endl;
        
        local_degrees_initialized = true;
    }
    
    /**
     * Get available degree for a node (from local copy)
     */
    int32_t get_local_available_degree(uint64_t node_id) const {
        initialize_local_degrees();
        auto it = local_degree_budgets.find(node_id);
        return (it != local_degree_budgets.end()) ? it->second : 0;
    }
    
    /**
     * Consume degree budget locally with incremental available list update
     */
    bool consume_local_degree(uint64_t node_id, int32_t amount = 1) const {
        initialize_local_degrees();
        auto it = local_degree_budgets.find(node_id);
        if (it != local_degree_budgets.end() && it->second >= amount) {
            it->second -= amount;
            
            // If budget drops to zero, remove from available list
            if (it->second == 0) {
                auto vec_it = std::find(local_available_nodes.begin(), 
                                       local_available_nodes.end(), node_id);
                if (vec_it != local_available_nodes.end()) {
                    local_available_nodes.erase(vec_it);
                }
            }
            
            return true;
        }
        return false;
    }
    
    /**
     * Get list of available nodes (maintained incrementally)
     */
    const std::vector<uint64_t>& get_local_available_nodes() const {
        initialize_local_degrees();
        return local_available_nodes;
    }
};

/**
 * Helper functions for stage implementations
 */
namespace degree_aware_stages {

/**
 * Get available nodes for this cluster (from local copy)
 */
std::unordered_set<uint64_t> get_cluster_available_nodes(
    const GraphTaskWithDegrees& task) {
    const auto& available_nodes = task.get_local_available_nodes();
    return std::unordered_set<uint64_t>(available_nodes.begin(), available_nodes.end());
}

/**
 * Select edge endpoint with preference for available degree nodes (using local copy)
 */
uint32_t select_edge_endpoint(
    const std::vector<uint32_t>& candidates,
    const Graph& g,
    const GraphTaskWithDegrees& task,
    std::mt19937& gen) {
    
    if (candidates.empty()) return UINT32_MAX;
    
    // Use local available nodes (no mutex needed!)
    const auto& available_nodes = task.get_local_available_nodes();
    std::unordered_set<uint64_t> available_set(available_nodes.begin(), available_nodes.end());
    
    // First pass: collect available candidates
    std::vector<uint32_t> available_candidates;
    for (uint32_t candidate : candidates) {
        uint64_t candidate_id = g.id_map[candidate];
        if (available_set.count(candidate_id)) {
            available_candidates.push_back(candidate);
        }
    }
    
    // Prefer available candidates
    if (!available_candidates.empty()) {
        std::uniform_int_distribution<size_t> dist(0, available_candidates.size() - 1);
        return available_candidates[dist(gen)];
    }
    
    // Fallback to any candidate
    std::uniform_int_distribution<size_t> dist(0, candidates.size() - 1);
    return candidates[dist(gen)];
}

} // namespace degree_aware_stages

#endif // AVAILABLE_NODE_DEGREES_H
