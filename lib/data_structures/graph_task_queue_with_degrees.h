#ifndef GRAPH_TASK_QUEUE_WITH_DEGREES_H
#define GRAPH_TASK_QUEUE_WITH_DEGREES_H

#include <nlohmann/json.hpp>
#include <queue>
#include <vector>
#include <memory>
#include <string>
#include <unordered_set>
#include <functional>
#include <mutex>
#include <atomic>
#include <omp.h>
#include "graph.h"
#include "clustering.h"
#include "available_node_degrees.h" // Our updated header with atomic manager
#include "../io/requirements_io.h"

using json = nlohmann::json;

class GraphTaskQueueWithDegrees {
private:
    std::queue<GraphTaskWithDegrees> task_queue;
    mutable std::mutex queue_mutex;
    std::atomic<size_t> tasks_processed{0};
    std::atomic<size_t> total_tasks_created{0};

    // Shared atomic degrees manager - used by ALL threads
    std::shared_ptr<AvailableNodeDegreesManager> degree_manager;
    
    // Store completed subgraphs for post-processing
    mutable std::mutex completed_subgraphs_mutex;
    std::vector<std::shared_ptr<Graph>> completed_subgraphs;
    
    // Task functions - now take the atomic task type
    std::function<void(GraphTaskWithDegrees&)> min_deg_enforce_fn;
    std::function<void(GraphTaskWithDegrees&)> cc_stitching_fn;
    std::function<void(GraphTaskWithDegrees&)> wcc_stitching_fn;

    // Map subgraph pointer to its mutex
    std::unordered_map<Graph*, std::unique_ptr<std::mutex>> subgraph_mutexes;
    std::mutex mutex_map_mutex;  // Protects the mutex map itself
    
    std::mutex& get_subgraph_mutex(Graph* subgraph) {
        std::lock_guard<std::mutex> lock(mutex_map_mutex);
        auto& mutex_ptr = subgraph_mutexes[subgraph];
        if (!mutex_ptr) {
            mutex_ptr = std::make_unique<std::mutex>();
        }
        return *mutex_ptr;
    }
    
    // Extract subgraph for a cluster (unchanged)
    std::shared_ptr<Graph> extract_subgraph(const Graph& original, 
                                           const std::unordered_set<uint32_t>& nodes, 
                                           const std::unordered_set<uint32_t>& missing_nodes = {}) {
        auto subgraph = std::make_shared<Graph>();
        
        uint32_t total_nodes = nodes.size() + missing_nodes.size();
        
        std::unordered_map<uint32_t, uint32_t> node_to_sub;
        std::unordered_map<uint64_t, uint32_t> missing_to_sub;
        uint32_t sub_idx = 0;
        
        for (uint32_t node : nodes) {
            node_to_sub[node] = sub_idx++;
        }
        
        for (uint64_t missing_node : missing_nodes) {
            missing_to_sub[missing_node] = sub_idx++;
        }
        
        subgraph->num_nodes = total_nodes;
        subgraph->row_ptr.resize(subgraph->num_nodes + 1, 0);
        
        std::vector<uint32_t> edge_counts(subgraph->num_nodes, 0);
        for (uint32_t orig_node : nodes) {
            uint32_t sub_node = node_to_sub[orig_node];
            
            for (uint32_t j = original.row_ptr[orig_node]; j < original.row_ptr[orig_node + 1]; j++) {
                uint32_t neighbor = original.col_idx[j];
                if (nodes.count(neighbor) > 0) {
                    edge_counts[sub_node]++;
                }
            }
        }
        
        for (uint32_t i = 0; i < subgraph->num_nodes; i++) {
            subgraph->row_ptr[i + 1] = subgraph->row_ptr[i] + edge_counts[i];
        }
        
        subgraph->col_idx.resize(subgraph->row_ptr.back());
        
        std::vector<uint32_t> current_pos(subgraph->num_nodes, 0);
        for (uint32_t i = 0; i < subgraph->num_nodes; i++) {
            current_pos[i] = subgraph->row_ptr[i];
        }
        
        for (uint32_t orig_node : nodes) {
            uint32_t sub_node = node_to_sub[orig_node];
            
            for (uint32_t j = original.row_ptr[orig_node]; j < original.row_ptr[orig_node + 1]; j++) {
                uint32_t neighbor = original.col_idx[j];
                if (nodes.count(neighbor) > 0) {
                    uint32_t sub_neighbor = node_to_sub[neighbor];
                    subgraph->col_idx[current_pos[sub_node]++] = sub_neighbor;
                }
            }
        }
        
        subgraph->node_map.clear();
        subgraph->id_map.resize(subgraph->num_nodes);
        
        for (const auto& [orig, sub] : node_to_sub) {
            if (original.id_map.size() > orig) {
                uint64_t original_id = original.id_map[orig];
                subgraph->node_map[original_id] = sub;
                subgraph->id_map[sub] = original_id;
            }
        }
        
        for (const auto& [missing_id, sub] : missing_to_sub) {
            subgraph->node_map[missing_id] = sub;
            subgraph->id_map[sub] = missing_id;
        }
        
        subgraph->num_edges = subgraph->col_idx.size() / 2;
        
        return subgraph;
    }

    bool try_get_task(GraphTaskWithDegrees& task) {
        std::lock_guard<std::mutex> lock(queue_mutex);
        if (task_queue.empty()) {
            return false;
        }
        
        task = task_queue.front();
        task_queue.pop();
        return true;
    }
    
    void add_task(const GraphTaskWithDegrees& task) {
        std::lock_guard<std::mutex> lock(queue_mutex);
        task_queue.push(task);
        total_tasks_created++;
    }
    
public:
    GraphTaskQueueWithDegrees() = default;
    
    /**
     * Initialize with degree deficits from JSON file - creates shared atomic manager
     */
    void initialize_degree_manager(const std::string& deficits_json_filename) {
        degree_manager = std::make_shared<AvailableNodeDegreesManager>(deficits_json_filename);
            
        auto stats = degree_manager->get_stats();
        std::cout << "Atomic degree manager initialized from JSON with " << stats.total_available_nodes 
                  << " available nodes, total budget: " << stats.total_available_degrees << std::endl;
    }
    
    void set_task_functions(
        std::function<void(GraphTaskWithDegrees&)> min_deg_enforce,
        std::function<void(GraphTaskWithDegrees&)> cc_stitch,
        std::function<void(GraphTaskWithDegrees&)> wcc_stitch) {
        
        min_deg_enforce_fn = min_deg_enforce;
        cc_stitching_fn = cc_stitch;
        wcc_stitching_fn = wcc_stitch;
    }
    
    void initialize_queue(const Graph& graph, const Clustering& clustering,
                         const ConnectivityRequirementsLoader& requirements) {
        if (!degree_manager) {
            throw std::runtime_error("Atomic degree manager must be initialized before queue initialization");
        }
        
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            while (!task_queue.empty()) {
                task_queue.pop();
            }
        }
        
        tasks_processed = 0;
        total_tasks_created = 0;
        
        // Create tasks for each cluster
        for (uint32_t cluster_idx = 0; cluster_idx < clustering.cluster_nodes.size(); cluster_idx++) {
            const auto& nodes = clustering.cluster_nodes[cluster_idx];
            const auto& missing_nodes = clustering.cluster_missing_nodes[cluster_idx];

            if (!missing_nodes.empty()) {
                std::cerr << "Warning: Cluster " << clustering.cluster_ids[cluster_idx]
                            << " has " << missing_nodes.size() 
                            << " missing nodes" << std::endl;
            }
            
            if (nodes.empty() && missing_nodes.empty()) {
                std::cerr << "Skipping empty cluster " << clustering.cluster_ids[cluster_idx] 
                          << " (index " << cluster_idx << ")" << std::endl;
                continue;
            }

            if (nodes.empty() && !missing_nodes.empty()) {
                std::cerr << "Warning: Cluster " << clustering.cluster_ids[cluster_idx] 
                          << " has no nodes but has missing nodes" << std::endl;
            }
            
            auto subgraph = extract_subgraph(graph, nodes, missing_nodes);
            subgraph->id = clustering.cluster_ids[cluster_idx];
            
            std::string cluster_id = clustering.cluster_ids[cluster_idx];
            uint32_t min_degree = requirements.get_connectivity_requirement(cluster_id);
            
            if (min_degree == 0) {
                std::cout << "Warning: No requirement found for cluster " << cluster_id 
                          << ", with " << nodes.size()
                          << " nodes, using default min_degree=1" << std::endl;
                min_degree = 1;
            }
            
            // Create cluster node IDs set for intersection operations
            std::unordered_set<uint64_t> cluster_node_ids;
            for (uint32_t node : nodes) {
                if (graph.id_map.size() > node) {
                    cluster_node_ids.insert(graph.id_map[node]);
                }
            }
            for (uint64_t missing_node : missing_nodes) {
                cluster_node_ids.insert(missing_node);
            }
            
            // Create task with shared atomic degree manager
            add_task(GraphTaskWithDegrees(
                subgraph,
                TaskType::MIN_DEG_ENFORCE,
                cluster_id,
                cluster_idx,
                min_degree,
                degree_manager,  // Shared across ALL tasks
                cluster_node_ids
            ));
        }
        
        std::cout << "Initialized task queue with " << task_queue.size() 
                  << " clusters with connectivity requirements" << std::endl;
    }
    
    bool process_next_task() {
        GraphTaskWithDegrees task(nullptr, TaskType::MIN_DEG_ENFORCE, "", 0, 0, nullptr, {});

        if (!try_get_task(task)) return false;

        int thread_id = omp_get_thread_num();
        
        // Lock this specific subgraph for the entire pipeline
        std::lock_guard<std::mutex> subgraph_lock(get_subgraph_mutex(task.subgraph.get()));
        
        std::cout << "Thread " << thread_id << " processing complete pipeline for cluster " 
                  << task.cluster_id << std::endl;

        // Process entire pipeline atomically with subgraph locked
        if (min_deg_enforce_fn) {
            min_deg_enforce_fn(task);
        }
        
        if (cc_stitching_fn) {
            cc_stitching_fn(task);
        }
        
        if (wcc_stitching_fn) {
            wcc_stitching_fn(task);
        }
        
        // Store completed subgraph (still under lock)
        {
            std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
            completed_subgraphs.push_back(task.subgraph);
        }

        tasks_processed++;
        return true;
        // subgraph_lock releases here
    }
    
    void process_all_tasks() {
        int num_threads = omp_get_max_threads();
        
        std::cout << "Starting work-stealing processing with " << num_threads 
                  << " threads using shared atomic degree manager" << std::endl;
        
        #pragma omp parallel num_threads(num_threads)
        {
            #pragma omp single
            {
                // Create tasks dynamically
                while (!is_empty()) {
                    #pragma omp task
                    {
                        process_next_task();
                    }
                }
            }
        }
        
        // Print final statistics
        auto final_stats = degree_manager->get_stats();
        std::cout << "Processed " << tasks_processed.load() << " tasks total" << std::endl;
        std::cout << "Final degree budget: " << final_stats.total_available_nodes 
                  << " nodes with " << final_stats.total_available_degrees 
                  << " remaining degrees" << std::endl;
    }
    
    // Access to atomic degree manager for external use
    std::shared_ptr<AvailableNodeDegreesManager> get_degree_manager() const {
        return degree_manager;
    }
    
    size_t queue_size() const {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return task_queue.size();
    }
    
    bool is_empty() const {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return task_queue.empty();
    }
    
    static std::string get_task_name(TaskType type) {
        switch (type) {
            case TaskType::MIN_DEG_ENFORCE: return "min_deg_enforce";
            case TaskType::CC_STITCHING: return "cc_stitching";
            case TaskType::WCC_STITCHING: return "wcc_stitching";
            default: return "unknown";
        }
    }

    void clear_completed_subgraphs() {
        std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
        completed_subgraphs.clear();
    }
    
    std::vector<std::shared_ptr<Graph>> get_completed_subgraphs() const {
        std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
        return completed_subgraphs;
    }
};

#endif // GRAPH_TASK_QUEUE_WITH_DEGREES_H
