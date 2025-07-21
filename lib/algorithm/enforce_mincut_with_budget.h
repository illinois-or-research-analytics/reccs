#ifndef ENFORCE_MINCUT_WITH_BUDGET_H
#define ENFORCE_MINCUT_WITH_BUDGET_H

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <random>
#include "../utils/statics.h"
#include "../data_structures/available_node_degrees.h"

/**
 * Degree-aware mincut enforcement using local degree management
 * Each task operates on its own local degree budgets - NO SHARED STATE!
 */
void enforce_mincut_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_cut_size = task.min_degree_requirement;

    auto start_time = std::chrono::steady_clock::now();
    const auto MAX_TIME = std::chrono::seconds(30); // 30 second timeout
    
    std::cout << "Starting local mincut enforcement on cluster "
              << g.id << ". Target minimum cut size: " << min_cut_size << std::endl;

    // Initialize local degrees for this task (lazy initialization)
    task.initialize_local_degrees();
    
    // Get local available nodes (no shared state!)
    const auto& local_available_nodes = task.get_local_available_nodes();
    
    std::cout << "Cluster has " << local_available_nodes.size() 
              << " nodes with local available degree budget" << std::endl;

    // Create lookup set for available nodes from LOCAL data
    std::unordered_set<uint64_t> available_node_set(local_available_nodes.begin(),
                                                   local_available_nodes.end());

    int iteration = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    // Pre-compute existing edges once for efficiency
    auto existing_edges = statics::compute_existing_edges(g);
    auto edge_exists = statics::create_edge_exists_checker(existing_edges);
    
    while (true) {
        iteration++;
        std::cout << "\n--- Iteration " << iteration << " ---" << std::endl;
        
        // Check timeout
        if (std::chrono::steady_clock::now() - start_time > MAX_TIME) {
            std::cout << "Timeout reached, stopping mincut enforcement" << std::endl;
            break;
        }
        
        // Get current minimum cut
        MincutResult result = compute_mincut(g);
        uint32_t current_cut_size = result.get_cut_size();

        std::cout << "[Cluster " << g.id << "] Current minimum cut size: " << current_cut_size << std::endl;

        // Check if we've reached the target
        if (current_cut_size >= min_cut_size) {
            std::cout << "[Cluster " << g.id << "] SUCCESS: Cut size " << current_cut_size 
                      << " meets the minimum requirement of " << min_cut_size 
                      << ". Stopping enforcement after " << iteration << " iterations." << std::endl;
            
            std::cout << "Total edges added: " << total_edges_added 
                      << ", Degree corrected: " << degree_corrected_edges 
                      << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
                      << "%)" << std::endl;
            return;
        }
        
        // Calculate how many edges we need to add
        uint32_t edges_needed = min_cut_size - current_cut_size;
        std::cout << "[Cluster " << g.id << "] Need to add " << edges_needed << " more edges between partitions." << std::endl;

        // Get the two partitions
        const std::vector<unsigned int>& light_partition = result.get_light_partition();
        const std::vector<unsigned int>& heavy_partition = result.get_heavy_partition();

        std::cout << "[Cluster " << g.id << "] Light partition size: " << light_partition.size() 
                  << ", Heavy partition size: " << heavy_partition.size() << std::endl;
        
        if (light_partition.empty() || heavy_partition.empty()) {
            std::cout << "[Cluster " << g.id << "] ERROR: One partition is empty. Cannot add cross-partition edges." << std::endl;
            return;
        }
        
        // OPTIMIZATION: Pre-compute all valid cross-partition edges
        std::cout << "Pre-filtering valid cross-partition edges..." << std::endl;
        
        std::vector<std::pair<uint32_t, uint32_t>> valid_budget_edges;  // Both nodes have budget
        std::vector<std::pair<uint32_t, uint32_t>> valid_mixed_edges;   // One node has budget  
        std::vector<std::pair<uint32_t, uint32_t>> valid_any_edges;     // Neither has budget
        
        for (unsigned int light_node : light_partition) {
            uint64_t light_node_id = g.id_map[light_node];
            bool light_has_budget = task.get_local_available_degree(light_node_id) > 0;
            
            for (unsigned int heavy_node : heavy_partition) {
                uint64_t heavy_node_id = g.id_map[heavy_node];
                bool heavy_has_budget = task.get_local_available_degree(heavy_node_id) > 0;
                
                uint32_t u = std::min(light_node, heavy_node);
                uint32_t v = std::max(light_node, heavy_node);
                
                // Skip if edge already exists
                if (edge_exists(u, v)) {
                    continue;
                }
                
                // Categorize by budget availability
                if (light_has_budget && heavy_has_budget) {
                    valid_budget_edges.emplace_back(u, v);
                } else if (light_has_budget || heavy_has_budget) {
                    valid_mixed_edges.emplace_back(u, v);
                } else {
                    valid_any_edges.emplace_back(u, v);
                }
            }
        }
        
        std::cout << "Found " << valid_budget_edges.size() << " edges with both nodes having budget, "
                  << valid_mixed_edges.size() << " edges with one node having budget, "
                  << valid_any_edges.size() << " edges with no budget" << std::endl;
        
        // Combine edge pools in priority order
        std::vector<std::pair<uint32_t, uint32_t>> candidate_edges;
        candidate_edges.reserve(valid_budget_edges.size() + valid_mixed_edges.size() + valid_any_edges.size());
        
        // Priority 1: Both nodes have budget
        candidate_edges.insert(candidate_edges.end(), valid_budget_edges.begin(), valid_budget_edges.end());
        
        // Priority 2: One node has budget
        candidate_edges.insert(candidate_edges.end(), valid_mixed_edges.begin(), valid_mixed_edges.end());
        
        // Priority 3: Neither node has budget (fallback)
        candidate_edges.insert(candidate_edges.end(), valid_any_edges.begin(), valid_any_edges.end());
        
        if (candidate_edges.empty()) {
            std::cout << "TERMINATION: No valid cross-partition edges available. "
                      << "All possible cross-partition edges already exist." << std::endl;
            std::cout << "Maximum achievable cut size with current partitioning: " << current_cut_size << std::endl;
            return;
        }
        
        // Randomly shuffle candidates to ensure randomness within priority groups
        if (!valid_budget_edges.empty()) {
            std::shuffle(candidate_edges.begin(), 
                        candidate_edges.begin() + valid_budget_edges.size(), gen);
        }
        if (!valid_mixed_edges.empty()) {
            std::shuffle(candidate_edges.begin() + valid_budget_edges.size(), 
                        candidate_edges.begin() + valid_budget_edges.size() + valid_mixed_edges.size(), gen);
        }
        if (!valid_any_edges.empty()) {
            std::shuffle(candidate_edges.begin() + valid_budget_edges.size() + valid_mixed_edges.size(), 
                        candidate_edges.end(), gen);
        }
        
        // Select edges to add (up to what we need and what's available)
        size_t edges_to_select = std::min(static_cast<size_t>(edges_needed), candidate_edges.size());
        std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
        
        std::cout << "Selecting " << edges_to_select << " edges from " << candidate_edges.size() 
                  << " candidates..." << std::endl;
        
        for (size_t i = 0; i < edges_to_select; ++i) {
            uint32_t u = candidate_edges[i].first;
            uint32_t v = candidate_edges[i].second;
            
            edges_to_add.emplace_back(u, v);
            total_edges_added++;
            
            // Track if this edge uses available degree budget using LOCAL data
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[v];
            if (task.get_local_available_degree(u_id) > 0 || task.get_local_available_degree(v_id) > 0) {
                degree_corrected_edges++;
            }
            
            // Consume local budgets immediately (no contention!)
            task.consume_local_degree(u_id, 1);
            task.consume_local_degree(v_id, 1);
            
            if ((i + 1) % 10 == 0 || i == edges_to_select - 1) {
                std::cout << "Selected cross-partition edge " << (i + 1) << "/" << edges_to_select 
                          << ": " << u << " <-> " << v << std::endl;
            }
        }
        
        std::cout << "Adding " << edges_to_add.size() << " new cross-partition edges..." << std::endl;
        
        // Batch add the edges to the graph
        add_edges_batch(g, edges_to_add);
        
        std::cout << "Added " << edges_to_add.size() << " edges in iteration " << iteration << std::endl;
        
        // Safety check to prevent infinite loops
        if (iteration > 100) {
            std::cout << "WARNING: Reached maximum iterations (100). "
                      << "Current cut size: " << current_cut_size 
                      << ", Target: " << min_cut_size << std::endl;
            break;
        }
    }
    
    // Report final local available nodes
    const auto& final_available = task.get_local_available_nodes();
    
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "Mincut enforcement completed in " << duration.count() << " seconds. "
              << "Total edges added: " << total_edges_added 
              << ", Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%). Remaining local budget nodes: " << final_available.size() << std::endl;
}

#endif // ENFORCE_MINCUT_WITH_BUDGET_H
