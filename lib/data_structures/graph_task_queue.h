#ifndef GRAPH_TASK_QUEUE_H
#define GRAPH_TASK_QUEUE_H

#include <queue>
#include <vector>
#include <memory>
#include <string>
#include <unordered_set>
#include <functional>
#include "graph.h"
#include "clustering.h"

enum class TaskType {
    DEGREE_ENFORCE = 0,
    CC_STITCHING = 1,
    WCC_STITCHING = 2,
    DEG_SEQ_MATCHING = 3
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
    
    // Constructor
    GraphTask(std::shared_ptr<Graph> g, TaskType t, const std::string& cid, uint32_t cidx)
        : subgraph(g), task_type(t), cluster_id(cid), cluster_idx(cidx) {}
};

class GraphTaskQueue {
private:
    // Verbose output flag
    bool verbose = false;

    // Queue to hold tasks
    std::queue<GraphTask> task_queue;
    
    // Task functions - to be implemented based on your algorithms
    std::function<void(Graph&)> degree_enforce_fn;
    std::function<void(Graph&)> cc_stitching_fn;
    std::function<void(Graph&)> wcc_stitching_fn;
    std::function<void(Graph&)> deg_seq_matching_fn;
    
    // Extract subgraph for a cluster
    std::shared_ptr<Graph> extract_subgraph(const Graph& original, 
                                           const std::unordered_set<uint32_t>& nodes) {
        auto subgraph = std::make_shared<Graph>();
        
        // Create node mapping from original to subgraph
        std::unordered_map<uint32_t, uint32_t> node_to_sub;
        uint32_t sub_idx = 0;
        for (uint32_t node : nodes) {
            node_to_sub[node] = sub_idx++;
        }
        
        subgraph->num_nodes = nodes.size();
        subgraph->row_ptr.resize(subgraph->num_nodes + 1, 0);
        
        // First pass: count edges for each node in subgraph
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
        
        // Build row_ptr
        for (uint32_t i = 0; i < subgraph->num_nodes; i++) {
            subgraph->row_ptr[i + 1] = subgraph->row_ptr[i] + edge_counts[i];
        }
        
        // Allocate col_idx
        subgraph->col_idx.resize(subgraph->row_ptr.back());
        
        // Second pass: fill col_idx
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
        
        // Store node mappings if needed
        subgraph->node_map.clear();
        subgraph->id_map.resize(subgraph->num_nodes);
        for (const auto& [orig, sub] : node_to_sub) {
            if (original.id_map.size() > orig) {
                uint64_t original_id = original.id_map[orig];
                subgraph->node_map[original_id] = sub;
                subgraph->id_map[sub] = original_id;
            }
        }
        
        // Count edges (each undirected edge is stored twice in CSR)
        subgraph->num_edges = subgraph->col_idx.size() / 2;
        
        return subgraph;
    }
    
public:
    // Constructor
    GraphTaskQueue() = default;

    void set_verbose(bool verbose) {
        this->verbose = verbose;
        if (verbose) {
            std::cout << "Verbose mode enabled for GraphTaskQueue" << std::endl;
        }
    }
    
    // Set task functions
    void set_task_functions(
        std::function<void(Graph&)> degree_enforce,
        std::function<void(Graph&)> cc_stitch,
        std::function<void(Graph&)> wcc_stitch,
        std::function<void(Graph&)> deg_seq_match) {
        
        degree_enforce_fn = degree_enforce;
        cc_stitching_fn = cc_stitch;
        wcc_stitching_fn = wcc_stitch;
        deg_seq_matching_fn = deg_seq_match;
    }
    
    // Initialize queue with all clusters
    void initialize_queue(const Graph& graph, const Clustering& clustering) {
        // Clear existing queue
        while (!task_queue.empty()) {
            task_queue.pop();
        }
        
        // Add each non-empty cluster to the queue
        for (uint32_t cluster_idx = 0; cluster_idx < clustering.cluster_nodes.size(); cluster_idx++) {
            const auto& nodes = clustering.cluster_nodes[cluster_idx];
            
            if (nodes.empty()) continue;
            
            // Extract subgraph for this cluster
            auto subgraph = extract_subgraph(graph, nodes);
            
            // Add to queue with first task
            task_queue.emplace(
                subgraph,
                TaskType::DEGREE_ENFORCE,
                clustering.cluster_ids[cluster_idx],
                cluster_idx
            );
        }
        
        if (verbose) {
            std::cout << "Initialized task queue with " << task_queue.size() 
                      << " tasks for clusters" << std::endl;
        }
    }
    
    // Process one task from the queue
    bool process_next_task() {
        if (task_queue.empty()) {
            return false;
        }
        
        GraphTask task = task_queue.front();
        task_queue.pop();
        
        // Process based on task type
        switch (task.task_type) {
            case TaskType::DEGREE_ENFORCE:
                if (verbose) {
                    std::cout << "Processing degree enforcement for cluster " << task.cluster_id << std::endl;
                }
                if (degree_enforce_fn) {
                    degree_enforce_fn(*task.subgraph);
                }
                // Add next task
                task.task_type = TaskType::CC_STITCHING;
                task_queue.push(task);
                break;
                
            case TaskType::CC_STITCHING:
                if (verbose) {
                    std::cout << "Processing CC stitching for cluster " << task.cluster_id << std::endl;
                }
                if (cc_stitching_fn) {
                    cc_stitching_fn(*task.subgraph);
                }
                // Add next task
                task.task_type = TaskType::WCC_STITCHING;
                task_queue.push(task);
                break;
                
            case TaskType::WCC_STITCHING:
                if (verbose) {
                    std::cout << "Processing WCC stitching for cluster " << task.cluster_id << std::endl;
                }
                if (wcc_stitching_fn) {
                    wcc_stitching_fn(*task.subgraph);
                }
                // Add next task
                task.task_type = TaskType::DEG_SEQ_MATCHING;
                task_queue.push(task);
                break;
                
            case TaskType::DEG_SEQ_MATCHING:
                if (verbose) {
                    std::cout << "Processing degree sequence matching for cluster " << task.cluster_id << std::endl;
                }
                if (deg_seq_matching_fn) {
                    deg_seq_matching_fn(*task.subgraph);
                }
                // Final task - don't re-queue
                if (verbose) {
                    std::cout << "Completed all tasks for cluster " << task.cluster_id << std::endl;
                }
                break;
        }
        
        return true;
    }
    
    // Process all tasks in the queue
    void process_all_tasks() {
        size_t tasks_processed = 0;
        while (process_next_task()) {
            tasks_processed++;
        }
        std::cout << "Processed " << tasks_processed << " tasks total" << std::endl;
    }
    
    // Get queue size
    size_t queue_size() const {
        return task_queue.size();
    }
    
    // Check if queue is empty
    bool is_empty() const {
        return task_queue.empty();
    }
    
    // Get task type name for debugging
    static std::string get_task_name(TaskType type) {
        switch (type) {
            case TaskType::DEGREE_ENFORCE: return "degree_enforce";
            case TaskType::CC_STITCHING: return "cc_stitching";
            case TaskType::WCC_STITCHING: return "wcc_stitching";
            case TaskType::DEG_SEQ_MATCHING: return "deg_seq_matching";
            default: return "unknown";
        }
    }
};

#endif // GRAPH_TASK_QUEUE_H
