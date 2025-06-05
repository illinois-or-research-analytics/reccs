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
#include "../io/requirements_io.h"

enum class TaskType {
    CONNECTIVITY_ENFORCE = 0,  // Combined degree + connectivity enforcement
    WCC_STITCHING = 1,
    DEG_SEQ_MATCHING = 2
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
    
    // Constructor
    GraphTask(std::shared_ptr<Graph> g, TaskType t, const std::string& cid, 
              uint32_t cidx, uint32_t min_deg = 1)
        : subgraph(g), task_type(t), cluster_id(cid), cluster_idx(cidx), 
          min_degree_requirement(min_deg) {}
};

class GraphTaskQueue {
private:
    std::queue<GraphTask> task_queue;
    
    // Task functions - to be implemented based on your algorithms
    std::function<void(Graph&, uint32_t)> connectivity_enforce_fn;
    std::function<void(Graph&)> wcc_stitching_fn;
    std::function<void(Graph&)> deg_seq_matching_fn;
    
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
    
public:
    // Constructor
    GraphTaskQueue() = default;
    
    // Set task functions
    void set_task_functions(
        std::function<void(Graph&, uint32_t)> connectivity_enforce,
        std::function<void(Graph&)> wcc_stitch,
        std::function<void(Graph&)> deg_seq_match) {
        
        connectivity_enforce_fn = connectivity_enforce;
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
            
            // Add to queue with first task (default min_degree = 1)
            task_queue.emplace(
                subgraph,
                TaskType::CONNECTIVITY_ENFORCE,
                clustering.cluster_ids[cluster_idx],
                cluster_idx,
                1  // default connectivity requirement
            );
        }
        
        std::cout << "Initialized task queue with " << task_queue.size() << " clusters" << std::endl;
    }
    
    // Initialize queue with all clusters and connectivity requirements
    void initialize_queue(const Graph& graph, const Clustering& clustering,
                         const ConnectivityRequirementsLoader& requirements) {
        // Clear existing queue
        while (!task_queue.empty()) {
            task_queue.pop();
        }
        
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
            
            if (nodes.empty()) continue;
            
            // Extract subgraph for this cluster
            auto subgraph = extract_subgraph(graph, nodes, missing_nodes);
            
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
            
            // Add to queue with first task
            task_queue.emplace(
                subgraph,
                TaskType::CONNECTIVITY_ENFORCE,
                cluster_id,
                cluster_idx,
                min_degree
            );
        }
        
        std::cout << "Initialized task queue with " << task_queue.size() 
                  << " clusters with connectivity requirements" << std::endl;
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
            case TaskType::CONNECTIVITY_ENFORCE:
                std::cout << "Processing connectivity enforcement for cluster " << task.cluster_id 
                          << " (min_degree=" << task.min_degree_requirement << ")" << std::endl;
                if (connectivity_enforce_fn) {
                    connectivity_enforce_fn(*task.subgraph, task.min_degree_requirement);
                }
                // Add next task
                task.task_type = TaskType::WCC_STITCHING;
                task_queue.push(task);
                break;
                
            case TaskType::WCC_STITCHING:
                std::cout << "Processing WCC stitching for cluster " << task.cluster_id << std::endl;
                if (wcc_stitching_fn) {
                    wcc_stitching_fn(*task.subgraph);
                }
                // Add next task
                task.task_type = TaskType::DEG_SEQ_MATCHING;
                task_queue.push(task);
                break;
                
            case TaskType::DEG_SEQ_MATCHING:
                std::cout << "Processing degree sequence matching for cluster " << task.cluster_id << std::endl;
                if (deg_seq_matching_fn) {
                    deg_seq_matching_fn(*task.subgraph);
                }
                // Final task - don't re-queue
                std::cout << "Completed all tasks for cluster " << task.cluster_id << std::endl;
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
            case TaskType::CONNECTIVITY_ENFORCE: return "connectivity_enforce";
            case TaskType::WCC_STITCHING: return "wcc_stitching";
            case TaskType::DEG_SEQ_MATCHING: return "deg_seq_matching";
            default: return "unknown";
        }
    }
};

#endif // GRAPH_TASK_QUEUE_H
