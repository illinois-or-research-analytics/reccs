#ifndef GRAPH_TASK_QUEUE_H
#define GRAPH_TASK_QUEUE_H

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
#include "../io/requirements_io.h"


using json = nlohmann::json;

enum class TaskType {
    MIN_DEG_ENFORCE = 0,
    CC_STITCHING = 1,
    WCC_STITCHING = 2
};

struct GraphTask {
    // The subgraph to process
    std::shared_ptr<Graph> subgraph;
    
    // Current task to perform
    TaskType task_type;
    
    // Cluster ID this subgraph belongs to
    std::string cluster_id;
    
    // Internal cluster index for efficiency
    uint32_t cluster_idx;
    
    // Minimum degree requirement for this cluster
    uint32_t min_degree_requirement;

    // Constructor - removed degree sequence parameters
    GraphTask(std::shared_ptr<Graph> g, TaskType t, const std::string& cid, 
              uint32_t cidx, uint32_t min_deg = 1)
        : subgraph(g), task_type(t), cluster_id(cid), cluster_idx(cidx), 
          min_degree_requirement(min_deg) {}
};

class GraphTaskQueue {
private:
    std::queue<GraphTask> task_queue;
    mutable std::mutex queue_mutex;  // Protects the queue
    std::atomic<size_t> tasks_processed{0};
    std::atomic<size_t> total_tasks_created{0};

    // Store completed subgraphs for post-processing
    mutable std::mutex completed_subgraphs_mutex;
    std::vector<std::shared_ptr<Graph>> completed_subgraphs;
    
    // Task functions - updated names to match TaskType enum
    std::function<void(Graph&, uint32_t)> min_deg_enforce_fn;
    std::function<void(Graph&, uint32_t)> cc_stitching_fn;
    std::function<void(Graph&, uint32_t)> wcc_stitching_fn;
    
    // Extract subgraph for a cluster
    std::shared_ptr<Graph> extract_subgraph(const Graph& original, 
                                           const std::unordered_set<uint32_t>& nodes, 
                                           const std::unordered_set<uint32_t>& missing_nodes = {}) {
        auto subgraph = std::make_shared<Graph>();
        
        // Total nodes = existing nodes + missing nodes
        uint32_t total_nodes = nodes.size() + missing_nodes.size();
        
        // Create node mapping from original to subgraph
        std::unordered_map<uint32_t, uint32_t> node_to_sub;
        std::unordered_map<uint64_t, uint32_t> missing_to_sub;
        uint32_t sub_idx = 0;
        
        // Map existing nodes first
        for (uint32_t node : nodes) {
            node_to_sub[node] = sub_idx++;
        }
        
        // Map missing nodes (they will be floaters)
        for (uint64_t missing_node : missing_nodes) {
            missing_to_sub[missing_node] = sub_idx++;
        }
        
        subgraph->num_nodes = total_nodes;
        subgraph->row_ptr.resize(subgraph->num_nodes + 1, 0);
        
        // First pass: count edges for each existing node in subgraph
        // Missing nodes will have 0 edges by default
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
        // Note: Missing nodes already have edge_counts[i] = 0, so they remain floaters
        
        // Build row_ptr
        for (uint32_t i = 0; i < subgraph->num_nodes; i++) {
            subgraph->row_ptr[i + 1] = subgraph->row_ptr[i] + edge_counts[i];
        }
        
        // Allocate col_idx
        subgraph->col_idx.resize(subgraph->row_ptr.back());
        
        // Second pass: fill col_idx for existing nodes only
        // Missing nodes have no edges to fill
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
        
        // Store node mappings for both existing and missing nodes
        subgraph->node_map.clear();
        subgraph->id_map.resize(subgraph->num_nodes);
        
        // Map existing nodes
        for (const auto& [orig, sub] : node_to_sub) {
            if (original.id_map.size() > orig) {
                uint64_t original_id = original.id_map[orig];
                subgraph->node_map[original_id] = sub;
                subgraph->id_map[sub] = original_id;
            }
        }
        
        // Map missing nodes (they use their original IDs directly)
        for (const auto& [missing_id, sub] : missing_to_sub) {
            subgraph->node_map[missing_id] = sub;
            subgraph->id_map[sub] = missing_id;
        }
        
        // Count edges (each undirected edge is stored twice in CSR)
        subgraph->num_edges = subgraph->col_idx.size() / 2;
        
        return subgraph;
    }

    // Thread-safe method to get next task
    bool try_get_task(GraphTask& task) {
        std::lock_guard<std::mutex> lock(queue_mutex);
        if (task_queue.empty()) {
            return false;
        }
        
        task = task_queue.front();
        task_queue.pop();
        return true;
    }
    
    // Thread-safe method to add task
    void add_task(const GraphTask& task) {
        std::lock_guard<std::mutex> lock(queue_mutex);
        task_queue.push(task);
        total_tasks_created++;
    }
    
public:
    // Constructor
    GraphTaskQueue() = default;
    
    // Set task functions - updated to match TaskType enum
    void set_task_functions(
        std::function<void(Graph&, uint32_t)> min_deg_enforce,
        std::function<void(Graph&, uint32_t)> cc_stitch,
        std::function<void(Graph&, uint32_t)> wcc_stitch) {
        
        min_deg_enforce_fn = min_deg_enforce;
        cc_stitching_fn = cc_stitch;
        wcc_stitching_fn = wcc_stitch;
    }
    
    // Initialize queue - removed degree_sequences_json parameter
    void initialize_queue(const Graph& graph, const Clustering& clustering,
                         const ConnectivityRequirementsLoader& requirements) {
        // Clear existing queue and reset counters
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            while (!task_queue.empty()) {
                task_queue.pop();
            }
        }
        
        tasks_processed = 0;
        total_tasks_created = 0;
        
        // Add each non-empty cluster to the queue
        for (uint32_t cluster_idx = 0; cluster_idx < clustering.cluster_nodes.size(); cluster_idx++) {
            const auto& nodes = clustering.cluster_nodes[cluster_idx];
            const auto& missing_nodes = clustering.cluster_missing_nodes[cluster_idx];

            // Check if cluster has any missing nodes
            if (!missing_nodes.empty()) {
                std::cerr << "Warning: Cluster " << clustering.cluster_ids[cluster_idx]
                            << " has " << missing_nodes.size() 
                            << " missing nodes" << std::endl;
            }
            
            if (nodes.empty() && missing_nodes.empty()) {
                std::cerr << "Skipping empty cluster " << clustering.cluster_ids[cluster_idx] 
                          << " (index " << cluster_idx << ")" << std::endl;
                continue;  // Skip empty clusters
            }

            if (nodes.empty() && !missing_nodes.empty()) {
                std::cerr << "Warning: Cluster " << clustering.cluster_ids[cluster_idx] 
                          << " has no nodes but has missing nodes" << std::endl;
            }
            
            // Extract subgraph for this cluster
            auto subgraph = extract_subgraph(graph, nodes, missing_nodes);

            // Assign cluster ID and index
            subgraph->id = clustering.cluster_ids[cluster_idx];
            
            // Get connectivity requirement for this cluster
            std::string cluster_id = clustering.cluster_ids[cluster_idx];
            uint32_t min_degree = requirements.get_connectivity_requirement(cluster_id);
            
            // Default to 1 if no requirement found
            if (min_degree == 0) {
                std::cout << "Warning: No requirement found for cluster " << cluster_id 
                          << ", with " << nodes.size()
                          << " nodes, using default min_degree=1" << std::endl;
                min_degree = 1;
            }
            
            // Add to queue with first task - updated to use MIN_DEG_ENFORCE
            add_task(GraphTask(
                subgraph,
                TaskType::MIN_DEG_ENFORCE,
                cluster_id,
                cluster_idx,
                min_degree
            ));
        }
        
        std::cout << "Initialized task queue with " << task_queue.size() 
                  << " clusters with connectivity requirements" << std::endl;
    }
    
    // Process one task from the queue
    bool process_next_task() {
        GraphTask task(nullptr, TaskType::MIN_DEG_ENFORCE, "", 0, 0);

        if (!try_get_task(task)) {
            return false;  // No more tasks
        }

        int thread_id = omp_get_thread_num();
        
        // Process based on task type - updated to match TaskType enum
        switch (task.task_type) {
            case TaskType::MIN_DEG_ENFORCE:
                std::cout << "Processing minimum degree enforcement for cluster " << task.cluster_id 
                          << " (min_degree=" << task.min_degree_requirement << ")" << std::endl;
                if (min_deg_enforce_fn) {
                    min_deg_enforce_fn(*task.subgraph, task.min_degree_requirement);
                }
                // Add next task
                task.task_type = TaskType::CC_STITCHING;
                add_task(task);
                break;
                
            case TaskType::CC_STITCHING:
                std::cout << "Processing CC stitching for cluster " << task.cluster_id << std::endl;
                if (cc_stitching_fn) {
                    cc_stitching_fn(*task.subgraph, task.min_degree_requirement);
                }
                // Add next task
                task.task_type = TaskType::WCC_STITCHING;
                add_task(task);
                break;
                
            case TaskType::WCC_STITCHING:
                std::cout << "Processing WCC stitching for cluster " << task.cluster_id << std::endl;
                if (wcc_stitching_fn) {
                    wcc_stitching_fn(*task.subgraph, task.min_degree_requirement);
                }
                
                // Store completed subgraph (final stage now)
                {
                    std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
                    completed_subgraphs.push_back(task.subgraph);
                }

                std::cout << "Completed all tasks for cluster " << task.cluster_id << std::endl;
                break;
        }
        // Increment processed tasks
        tasks_processed++;
        
        return true;
    }
    
    // Process all tasks in the queue
    void process_all_tasks() {
        int num_threads = omp_get_max_threads();
        
        std::cout << "Starting work-stealing processing with " << num_threads << " threads" << std::endl;
        
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
        
        std::cout << "Processed " << tasks_processed.load() << " tasks total" << std::endl;
    }
    
    // Get queue size
    size_t queue_size() const {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return task_queue.size();
    }
    
    // Check if queue is empty
    bool is_empty() const {
        std::lock_guard<std::mutex> lock(queue_mutex);
        return task_queue.empty();
    }
    
    // Get task type name for debugging - updated to match TaskType enum
    static std::string get_task_name(TaskType type) {
        switch (type) {
            case TaskType::MIN_DEG_ENFORCE: return "min_deg_enforce";
            case TaskType::CC_STITCHING: return "cc_stitching";
            case TaskType::WCC_STITCHING: return "wcc_stitching";
            default: return "unknown";
        }
    }

    // Clear completed subgraphs storage
    void clear_completed_subgraphs() {
        std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
        completed_subgraphs.clear();
    }
    
    // Get access to completed subgraphs for post-processing
    std::vector<std::shared_ptr<Graph>> get_completed_subgraphs() const {
        std::lock_guard<std::mutex> lock(completed_subgraphs_mutex);
        return completed_subgraphs;
    }
};

#endif // GRAPH_TASK_QUEUE_H
