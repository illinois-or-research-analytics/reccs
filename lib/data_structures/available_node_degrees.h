#ifndef AVAILABLE_NODE_DEGREES_H
#define AVAILABLE_NODE_DEGREES_H

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <fstream>
#include <iostream>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <vector>
#include "graph.h"

using json = nlohmann::json;

/**
 * Thread-safe atomic degree manager - shared across all threads
 * Uses atomic operations for degree consumption to prevent race conditions
 */
class AvailableNodeDegreesManager {
private:
    // Thread-safe degree budgets using atomic integers
    // Map structure is immutable after initialization, so no mutex needed for map access
    std::unordered_map<uint64_t, std::atomic<int32_t>> global_degree_budgets;
    
    // Cache of available node IDs (updated periodically during runtime)
    mutable std::vector<uint64_t> cached_available_nodes;
    mutable std::mutex cache_mutex;
    mutable std::atomic<bool> cache_dirty{true};

public:
    /**
     * Constructor that loads degree deficits from JSON file
     */
    AvailableNodeDegreesManager(const std::string& json_filename) {
        load_degree_budgets_from_json(json_filename);
    }
    
    /**
     * Load degree budgets from JSON file - called once during initialization
     * NO SYNCHRONIZATION NEEDED - called before threading starts
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
        // No synchronization needed - we're single-threaded during initialization
        for (const auto& [node_id_str, deficit] : j.items()) {
            try {
                uint64_t node_id = std::stoull(node_id_str);
                int32_t budget = deficit.get<int32_t>();
                
                // Only store positive budgets (matching original Python logic)
                if (budget > 0) {
                    global_degree_budgets.emplace(
                        std::piecewise_construct,
                        std::forward_as_tuple(node_id),
                        std::forward_as_tuple(budget)
                    );
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
        
        // Initialize cache - no synchronization needed during init
        update_available_nodes_cache_unsafe();
    }
    
    /**
     * Thread-safe degree consumption with atomic operations
     * Returns true if consumption was successful, false if insufficient budget
     */
    bool consume_degree(uint64_t node_id, int32_t amount = 1) {
        // No mutex needed for map access - map structure is immutable after init
        auto it = global_degree_budgets.find(node_id);
        if (it == global_degree_budgets.end()) {
            return false; // Node not found or has no budget
        }
        
        // Atomic compare-and-swap loop to safely consume degree
        int32_t current_budget = it->second.load();
        while (current_budget >= amount) {
            if (it->second.compare_exchange_weak(current_budget, current_budget - amount)) {
                // Successfully consumed the degree
                cache_dirty.store(true); // Mark cache as needing update
                return true;
            }
            // compare_exchange_weak failed, current_budget now contains the actual current value
            // Loop will retry with the updated current_budget
        }
        
        return false; // Insufficient budget
    }
    
    /**
     * Get current available degree for a node (thread-safe read)
     */
    int32_t get_available_degree(uint64_t node_id) const {
        // No mutex needed - map structure is immutable after init
        auto it = global_degree_budgets.find(node_id);
        return (it != global_degree_budgets.end()) ? it->second.load() : 0;
    }
    
    /**
     * Get list of nodes that currently have positive budget
     * Uses cached version for performance, updates cache when dirty
     */
    std::vector<uint64_t> get_available_nodes() const {
        // Check if cache needs updating
        if (cache_dirty.load()) {
            update_available_nodes_cache();
        }
        
        std::lock_guard<std::mutex> cache_lock(cache_mutex);
        return cached_available_nodes;
    }
    
    /**
     * Get available nodes for a specific cluster (intersection with cluster nodes)
     */
    std::vector<uint64_t> get_available_nodes_for_cluster(
        const std::unordered_set<uint64_t>& cluster_node_ids) const {
        
        std::vector<uint64_t> available_nodes;
        
        // No mutex needed for map iteration - map structure is immutable after init
        for (uint64_t node_id : cluster_node_ids) {
            auto it = global_degree_budgets.find(node_id);
            if (it != global_degree_budgets.end() && it->second.load() > 0) {
                available_nodes.push_back(node_id);
            }
        }
        
        return available_nodes;
    }
    
    /**
     * Try to consume degree with timeout/retry logic for high contention scenarios
     */
    bool consume_degree_with_retry(uint64_t node_id, int32_t amount = 1, int max_retries = 100) {
        for (int retry = 0; retry < max_retries; ++retry) {
            if (consume_degree(node_id, amount)) {
                return true;
            }
            
            // Small backoff on high contention
            if (retry > 10) {
                std::this_thread::yield();
            }
        }
        return false;
    }
    
    /**
     * Get statistics about available degrees (thread-safe)
     */
    struct AvailableStats {
        size_t total_available_nodes;
        int32_t total_available_degrees;
        double avg_available_degree;
    };
    
    AvailableStats get_stats() const {
        // No mutex needed for map iteration - map structure is immutable after init
        AvailableStats stats;
        stats.total_available_nodes = 0;
        stats.total_available_degrees = 0;
        
        for (const auto& [node_id, atomic_budget] : global_degree_budgets) {
            int32_t current_budget = atomic_budget.load();
            if (current_budget > 0) {
                stats.total_available_nodes++;
                stats.total_available_degrees += current_budget;
            }
        }
        
        stats.avg_available_degree = (stats.total_available_nodes > 0) ? 
            static_cast<double>(stats.total_available_degrees) / stats.total_available_nodes : 0.0;
            
        return stats;
    }

private:
    /**
     * Update the cached list of available nodes (thread-safe version)
     */
    void update_available_nodes_cache() const {
        std::vector<uint64_t> new_available_nodes;
        
        // No mutex needed for map iteration - map structure is immutable after init
        for (const auto& [node_id, atomic_budget] : global_degree_budgets) {
            if (atomic_budget.load() > 0) {
                new_available_nodes.push_back(node_id);
            }
        }
        
        {
            std::lock_guard<std::mutex> cache_lock(cache_mutex);
            cached_available_nodes = std::move(new_available_nodes);
        }
        
        cache_dirty.store(false);
    }
    
    /**
     * Update cache during initialization (no synchronization needed)
     */
    void update_available_nodes_cache_unsafe() {
        cached_available_nodes.clear();
        
        for (const auto& [node_id, atomic_budget] : global_degree_budgets) {
            if (atomic_budget.load() > 0) {
                cached_available_nodes.push_back(node_id);
            }
        }
        
        cache_dirty.store(false);
    }
};

/**
 * Simplified GraphTask - no more local degree management needed!
 */
struct GraphTaskWithDegrees {
    // Original task data
    std::shared_ptr<Graph> subgraph;
    TaskType task_type;
    std::string cluster_id;
    uint32_t cluster_idx;
    uint32_t min_degree_requirement;
    
    // Reference to shared atomic degree manager
    std::shared_ptr<AvailableNodeDegreesManager> degree_manager;
    
    // Cluster node IDs
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
    
    /**
     * Thread-safe degree consumption - directly uses shared manager
     */
    bool consume_local_degree(uint64_t node_id, int32_t amount = 1) const {
        return degree_manager->consume_degree(node_id, amount);
    }
    
    /**
     * Get available degree for a node
     */
    int32_t get_local_available_degree(uint64_t node_id) const {
        return degree_manager->get_available_degree(node_id);
    }
    
    /**
     * Get list of available nodes for this cluster
     */
    std::vector<uint64_t> get_local_available_nodes() const {
        return degree_manager->get_available_nodes_for_cluster(cluster_node_ids);
    }
};

/**
 * Updated helper functions for stage implementations
 */
namespace degree_aware_stages {

/**
 * Get available nodes for this cluster (from shared atomic manager)
 */
std::unordered_set<uint64_t> get_cluster_available_nodes(
    const GraphTaskWithDegrees& task) {
    const auto available_nodes = task.get_local_available_nodes();
    return std::unordered_set<uint64_t>(available_nodes.begin(), available_nodes.end());
}

/**
 * Select edge endpoint with preference for available degree nodes (thread-safe)
 */
uint32_t select_edge_endpoint(
    const std::vector<uint32_t>& candidates,
    const Graph& g,
    const GraphTaskWithDegrees& task,
    std::mt19937& gen) {
    
    if (candidates.empty()) return UINT32_MAX;
    
    // Get current available nodes from shared manager
    const auto available_nodes_vec = task.get_local_available_nodes();
    std::unordered_set<uint64_t> available_set(available_nodes_vec.begin(), available_nodes_vec.end());
    
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
