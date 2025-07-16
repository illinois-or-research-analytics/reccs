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
    // Dense storage for available degrees (node_id -> available_budget)
    std::unique_ptr<std::atomic<int32_t>[]> available_degrees;  // FIXED: Use unique_ptr array
    size_t available_degrees_size;  // FIXED: Track size separately

    // Reference degree sequence (empirical target degrees)
    std::shared_ptr<const std::vector<uint32_t>> reference_degrees;

    // Pre-computed available nodes per cluster to avoid repeated intersections
    std::unordered_map<std::string, std::vector<uint64_t>> cluster_available_cache;
    mutable std::mutex cache_mutex;
    
    // Reference graphs for degree computation
    const Graph* synthetic_graph;

public:
    AvailableNodeDegreesManager(
        const Graph& syn_graph,
        std::shared_ptr<const std::vector<uint32_t>> ref_degrees)
        : reference_degrees(ref_degrees), synthetic_graph(&syn_graph),
          available_degrees_size(0) {  // Initialize size to 0
        
        initialize_available_degrees();
    }
    
    /**
     * Initialize available degrees based on empirical vs synthetic degree differences
     * Replicates: deg_diff = emp_degrees[c_node] - sbm_degrees[c_node]
     */
    void initialize_available_degrees() {
        size_t max_nodes = std::max(reference_degrees->size(), synthetic_graph->num_nodes);
        
        // Allocate array of atomics using unique_ptr
        available_degrees_size = max_nodes;
        available_degrees = std::make_unique<std::atomic<int32_t>[]>(max_nodes);  // FIXED: Now matches declaration
        
        // Initialize all to zero
        for (size_t i = 0; i < max_nodes; ++i) {
            available_degrees[i].store(0);
        }

        size_t positive_degrees = 0;
        
        // For each node in the synthetic graph, compute available budget
        for (uint32_t i = 0; i < synthetic_graph->num_nodes; ++i) {
            // Get the original node ID
            uint64_t original_id = synthetic_graph->id_map[i];
            
            // Use original_id as index into reference degrees (assuming dense mapping)
            if (original_id < reference_degrees->size()) {
                uint32_t reference_degree = (*reference_degrees)[original_id];
                uint32_t current_degree = synthetic_graph->get_degree(i);
                
                int32_t budget = static_cast<int32_t>(reference_degree) - static_cast<int32_t>(current_degree);
                
                if (budget > 0 && original_id < available_degrees_size) {  // FIXED: Use size variable
                    available_degrees[original_id].store(budget);
                    positive_degrees++;
                }
            }
        }
        
        std::cout << "Initialized " << positive_degrees 
                  << " nodes with available degree budget from reference sequence" << std::endl;
    }

    /**
     * Pre-compute available nodes for each cluster
     */
    void precompute_cluster_available_nodes(
        const std::string& cluster_id,
        const std::unordered_set<uint64_t>& cluster_nodes) {
        
        std::vector<uint64_t> available_nodes;
        available_nodes.reserve(cluster_nodes.size());
        
        for (uint64_t node_id : cluster_nodes) {
            if (get_available_degree(node_id) > 0) {
                available_nodes.push_back(node_id);
            }
        }
        
        std::lock_guard<std::mutex> lock(cache_mutex);
        cluster_available_cache[cluster_id] = std::move(available_nodes);
    }

    /**
     * Get pre-computed available nodes for a cluster
     */
    std::vector<uint64_t> get_cluster_available_nodes(const std::string& cluster_id) const {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cluster_available_cache.find(cluster_id);
        return (it != cluster_available_cache.end()) ? it->second : std::vector<uint64_t>{};
    }
    
    /**
     * Fast degree lookup - O(1) access
     */
    int32_t get_available_degree(uint64_t node_id) const {
        if (node_id < available_degrees_size) {  // FIXED: Use size variable
            return available_degrees[node_id].load(std::memory_order_relaxed);
        }
        return 0;
    }
    
    /**
     * Try to consume degree budget
     */
    bool try_consume_degree(uint64_t node_id, int32_t amount = 1) {
        if (node_id >= available_degrees_size) return false;  // FIXED: Use size variable
        
        int32_t current = available_degrees[node_id].load(std::memory_order_relaxed);
        while (current >= amount) {
            if (available_degrees[node_id].compare_exchange_weak(
                current, current - amount, std::memory_order_relaxed)) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Batch consume multiple degrees
     */
    bool try_consume_degrees_batch(const std::vector<std::pair<uint64_t, int32_t>>& consumptions) {
        // First pass: check if all consumptions are possible
        for (const auto& [node_id, amount] : consumptions) {
            if (get_available_degree(node_id) < amount) {
                return false;
            }
        }
        
        // Second pass: perform all consumptions
        for (const auto& [node_id, amount] : consumptions) {
            if (!try_consume_degree(node_id, amount)) {
                // Race condition - some other thread consumed budget
                // This is rare, but we handle it gracefully
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * Update cluster cache after degree changes
     */
    void update_cluster_cache(const std::string& cluster_id, 
                             const std::unordered_set<uint64_t>& cluster_nodes) {
        std::vector<uint64_t> updated_available;
        
        {
            std::lock_guard<std::mutex> lock(cache_mutex);
            auto it = cluster_available_cache.find(cluster_id);
            if (it == cluster_available_cache.end()) return;
            
            // Filter out nodes with zero available degrees
            for (uint64_t node_id : it->second) {
                if (get_available_degree(node_id) > 0) {
                    updated_available.push_back(node_id);
                }
            }
            
            it->second = std::move(updated_available);
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
        AvailableStats stats;
        stats.total_available_nodes = 0;
        stats.total_available_degrees = 0;
        
        // FIXED: Use size variable and array indexing
        for (size_t i = 0; i < available_degrees_size; ++i) {
            int32_t val = available_degrees[i].load(std::memory_order_relaxed);
            if (val > 0) {
                stats.total_available_nodes++;
                stats.total_available_degrees += val;
            }
        }
        
        stats.avg_available_degree = (stats.total_available_nodes > 0) ? 
            static_cast<double>(stats.total_available_degrees) / stats.total_available_nodes : 0.0;
            
        return stats;
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
    return std::unordered_set<uint64_t>(task.degree_manager->get_cluster_available_nodes(task.cluster_id).begin(),
                                       task.degree_manager->get_cluster_available_nodes(task.cluster_id).end());
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
    
    // Get cached available nodes
    auto cluster_available_nodes = task.degree_manager->get_cluster_available_nodes(task.cluster_id);
    std::unordered_set<uint64_t> available_set(cluster_available_nodes.begin(), 
                                              cluster_available_nodes.end());
    
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
