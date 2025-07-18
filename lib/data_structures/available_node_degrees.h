#ifndef AVAILABLE_NODE_DEGREES_H
#define AVAILABLE_NODE_DEGREES_H

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "graph.h"

/**
 * Global immutable degree mapping - computed once, never modified
 * Each task gets its own local copy of relevant degrees
 */
class AvailableNodeDegreesManager {
private:
    // Global immutable mapping: node_id -> available_degree_budget
    std::unordered_map<uint64_t, int32_t> global_degree_budgets;
    
    // Reference degree sequence (for initialization)
    std::shared_ptr<const std::vector<uint32_t>> reference_degrees;
    
    // Reference to synthetic graph (for initialization)
    const Graph* synthetic_graph;

public:
    AvailableNodeDegreesManager(
        const Graph& syn_graph,
        std::shared_ptr<const std::vector<uint32_t>> ref_degrees)
        : reference_degrees(ref_degrees), synthetic_graph(&syn_graph) {
        
        initialize_global_degree_budgets();
    }
    
    /**
     * Initialize global degree budgets once - never modified after this
     */
    void initialize_global_degree_budgets() {
        std::cout << "Computing global degree budgets..." << std::endl;
        
        // Collect all nodes and their synthetic degrees
        std::unordered_map<uint64_t, uint32_t> synthetic_degrees;
        for (uint32_t i = 0; i < synthetic_graph->num_nodes; ++i) {
            uint64_t node_id = synthetic_graph->id_map[i];
            synthetic_degrees[node_id] = synthetic_graph->get_degree(i);
        }
        
        // Add missing nodes (not in synthetic graph) as degree 0
        for (size_t i = 0; i < reference_degrees->size(); ++i) {
            uint64_t node_id = i;  // Assuming dense node IDs
            if (synthetic_degrees.find(node_id) == synthetic_degrees.end()) {
                synthetic_degrees[node_id] = 0;  // Missing node has degree 0
            }
        }
        
        // Sort ALL nodes by synthetic degree (descending) for proper sequence matching
        std::vector<std::pair<uint32_t, uint64_t>> sorted_all_degrees;  // (degree, node_id)
        for (const auto& [node_id, degree] : synthetic_degrees) {
            sorted_all_degrees.emplace_back(degree, node_id);
        }
        
        std::sort(sorted_all_degrees.begin(), sorted_all_degrees.end(), 
                  [](const auto& a, const auto& b) {
                      return a.first > b.first;  // Sort by degree descending
                  });
        
        // Match with reference sequence and compute global budgets
        size_t positive_budgets = 0;
        for (size_t i = 0; i < sorted_all_degrees.size() && i < reference_degrees->size(); ++i) {
            uint64_t node_id = sorted_all_degrees[i].second;
            uint32_t synthetic_degree = sorted_all_degrees[i].first;
            uint32_t reference_degree = (*reference_degrees)[i];
            
            int32_t budget = static_cast<int32_t>(reference_degree) - static_cast<int32_t>(synthetic_degree);
            
            // Store budget (can be negative, zero, or positive)
            global_degree_budgets[node_id] = budget;
            
            if (budget > 0) {
                positive_budgets++;
            }
            
            // Debug output for first few nodes
            if (i < 10) {
                std::cout << "  Node " << node_id << ": ref=" << reference_degree 
                          << ", syn=" << synthetic_degree << ", budget=" << budget << std::endl;
            }
        }
        
        std::cout << "Global degree budgets computed: " << global_degree_budgets.size() 
                  << " total nodes, " << positive_budgets << " with positive budget" << std::endl;
    }
    
    /**
     * Create local degree budget map for a specific cluster
     * Each task gets its own independent copy - NO SHARED STATE!
     */
    std::unordered_map<uint64_t, int32_t> create_local_degrees(
        const std::unordered_set<uint64_t>& cluster_node_ids) const {
        
        std::unordered_map<uint64_t, int32_t> local_budgets;
        
        // Extract only the nodes belonging to this cluster
        for (uint64_t node_id : cluster_node_ids) {
            auto it = global_degree_budgets.find(node_id);
            if (it != global_degree_budgets.end()) {
                local_budgets[node_id] = it->second;  // Copy the budget
            } else {
                local_budgets[node_id] = 0;  // Default to 0 if not found
            }
        }
        
        std::cout << "Created local degrees for cluster with " << local_budgets.size() 
                  << " nodes" << std::endl;
        
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
            if (it != global_degree_budgets.end() && it->second > 0) {
                available_nodes.push_back(node_id);
            }
        }
        
        return available_nodes;
    }
    
    /**
     * Legacy methods for compatibility (not needed with local degrees approach)
     */
    void precompute_cluster_available_nodes(
        const std::string& cluster_id,
        const std::unordered_set<uint64_t>& cluster_nodes) {
        // No-op - each task creates its own local degrees
    }
    
    std::vector<uint64_t> get_cluster_available_nodes(const std::string& cluster_id) const {
        // No-op - each task has its own local degrees
        return {};
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
        stats.total_available_nodes = 0;
        stats.total_available_degrees = 0;
        
        for (const auto& [node_id, budget] : global_degree_budgets) {
            if (budget > 0) {
                stats.total_available_nodes++;
                stats.total_available_degrees += budget;
            }
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
