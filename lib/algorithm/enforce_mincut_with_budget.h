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
 * Degree-aware mincut enforcement
 * Replicates Python logic for Stage 3: ensuring minimum cut while 
 * prioritizing nodes with available degree budget
 */
void enforce_mincut_with_budget(GraphTaskWithDegrees& task) {
    Graph& g = *task.subgraph;
    uint32_t min_cut_size = task.min_degree_requirement;
    
    std::cout << "Starting mincut enforcement with budget on cluster "
              << g.id << ". Target minimum cut size: " << min_cut_size << std::endl;

    // Get available nodes for this cluster
    auto cluster_available_nodes = task.degree_manager->get_cluster_available_nodes(task.cluster_id);
    
    std::cout << "Cluster has " << cluster_available_nodes.size() 
              << " nodes with available degree budget" << std::endl;

    // Create lookup set for available nodes
    std::unordered_set<uint64_t> available_node_set(cluster_available_nodes.begin(),
                                                   cluster_available_nodes.end());

    int iteration = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Counters for statistics
    size_t total_edges_added = 0;
    size_t degree_corrected_edges = 0;
    
    while (true) {
        iteration++;
        std::cout << "\n--- Iteration " << iteration << " ---" << std::endl;
        
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
        
        // Track edges to add (avoid duplicates and existing edges)
        std::vector<std::pair<uint32_t, uint32_t>> edges_to_add;
        
        // Use statics.h for efficient edge existence checking
        auto existing_edges = statics::compute_existing_edges(g);
        auto edge_exists = statics::create_edge_exists_checker(existing_edges);
        
        // Categorize nodes by available budget (Python-style)
        std::vector<uint32_t> light_available, light_all;
        std::vector<uint32_t> heavy_available, heavy_all;
        
        // Build light partition candidate lists
        for (unsigned int node : light_partition) {
            light_all.push_back(node);
            
            uint64_t node_id = g.id_map[node];
            if (available_node_set.count(node_id) && 
                task.degree_manager->get_available_degree(node_id) > 0) {
                light_available.push_back(node);
            }
        }
        
        // Build heavy partition candidate lists
        for (unsigned int node : heavy_partition) {
            heavy_all.push_back(node);
            
            uint64_t node_id = g.id_map[node];
            if (available_node_set.count(node_id) && 
                task.degree_manager->get_available_degree(node_id) > 0) {
                heavy_available.push_back(node);
            }
        }

        std::cout << "Light partition: " << light_available.size() 
                  << " with budget, " << light_all.size() - light_available.size() << " without" << std::endl;
        std::cout << "Heavy partition: " << heavy_available.size() 
                  << " with budget, " << heavy_all.size() - heavy_available.size() << " without" << std::endl;
        
        // Add edges using Python-style random selection
        for (uint32_t edge_count = 0; edge_count < edges_needed; ++edge_count) {
            uint32_t selected_light = UINT32_MAX;
            uint32_t selected_heavy = UINT32_MAX;
            bool found_edge = false;
            
            // Python-style selection priority:
            // Priority 1: Both sides have available budget
            if (!light_available.empty() && !heavy_available.empty()) {
                std::uniform_int_distribution<size_t> dist_light(0, light_available.size() - 1);
                std::uniform_int_distribution<size_t> dist_heavy(0, heavy_available.size() - 1);
                selected_light = light_available[dist_light(gen)];
                selected_heavy = heavy_available[dist_heavy(gen)];
                found_edge = true;
            }
            // Priority 2: Light has budget, heavy any
            else if (!light_available.empty() && !heavy_all.empty()) {
                std::uniform_int_distribution<size_t> dist_light(0, light_available.size() - 1);
                std::uniform_int_distribution<size_t> dist_heavy(0, heavy_all.size() - 1);
                selected_light = light_available[dist_light(gen)];
                selected_heavy = heavy_all[dist_heavy(gen)];
                found_edge = true;
            }
            // Priority 2: Light any, heavy has budget
            else if (!light_all.empty() && !heavy_available.empty()) {
                std::uniform_int_distribution<size_t> dist_light(0, light_all.size() - 1);
                std::uniform_int_distribution<size_t> dist_heavy(0, heavy_available.size() - 1);
                selected_light = light_all[dist_light(gen)];
                selected_heavy = heavy_available[dist_heavy(gen)];
                found_edge = true;
            }
            // Priority 3: Neither side has budget (fallback)
            else if (!light_all.empty() && !heavy_all.empty()) {
                std::uniform_int_distribution<size_t> dist_light(0, light_all.size() - 1);
                std::uniform_int_distribution<size_t> dist_heavy(0, heavy_all.size() - 1);
                selected_light = light_all[dist_light(gen)];
                selected_heavy = heavy_all[dist_heavy(gen)];
                found_edge = true;
            }
            
            if (!found_edge) {
                std::cout << "Warning: Could not find nodes for cross-partition edge " 
                          << (edge_count + 1) << "/" << edges_needed << std::endl;
                break;
            }
            
            // Check if edge already exists
            uint32_t u = std::min(selected_light, selected_heavy);
            uint32_t v = std::max(selected_light, selected_heavy);
            
            if (edge_exists(u, v)) {
                // Skip existing edge, but don't count against our limit
                edge_count--;
                continue;
            }
            
            // Check if already in our planned edges
            bool already_planned = false;
            for (const auto& edge : edges_to_add) {
                if (edge.first == u && edge.second == v) {
                    already_planned = true;
                    break;
                }
            }
            
            if (already_planned) {
                // Skip already planned edge, but don't count against our limit
                edge_count--;
                continue;
            }
            
            // Add the edge
            edges_to_add.emplace_back(u, v);
            total_edges_added++;
            
            // Track if this edge uses available degree budget
            uint64_t u_id = g.id_map[u];
            uint64_t v_id = g.id_map[v];
            if (available_node_set.count(u_id) || available_node_set.count(v_id)) {
                degree_corrected_edges++;
            }
            
            std::cout << "Cross-partition edge " << (edge_count + 1) << "/" << edges_needed 
                      << ": " << selected_light << " <-> " << selected_heavy << std::endl;
        }
        
        if (edges_to_add.empty()) {
            std::cout << "TERMINATION: No more edges can be added between partitions. "
                      << "All possible cross-partition edges already exist." << std::endl;
            std::cout << "Maximum achievable cut size with current partitioning: " << current_cut_size << std::endl;
            return;
        }
        
        std::cout << "Adding " << edges_to_add.size() << " new cross-partition edges..." << std::endl;
        
        // Batch add the edges to the graph
        add_edges_batch(g, edges_to_add);
        
        // Batch consume degree budgets efficiently
        std::vector<std::pair<uint64_t, int32_t>> consumptions;
        consumptions.reserve(edges_to_add.size() * 2);
        
        for (const auto& edge : edges_to_add) {
            uint64_t u_id = g.id_map[edge.first];
            uint64_t v_id = g.id_map[edge.second];
            
            if (available_node_set.count(u_id)) {
                consumptions.emplace_back(u_id, 1);
            }
            if (available_node_set.count(v_id)) {
                consumptions.emplace_back(v_id, 1);
            }
        }
        
        // Try batch consumption
        task.degree_manager->try_consume_degrees_batch(consumptions);
        
        // Update cluster cache
        task.degree_manager->update_cluster_cache(task.cluster_id, task.cluster_node_ids);
        
        std::cout << "Added " << edges_to_add.size() << " edges in iteration " << iteration << std::endl;
        
        // Safety check to prevent infinite loops
        if (iteration > 100) {
            std::cout << "WARNING: Reached maximum iterations (100). "
                      << "Current cut size: " << current_cut_size 
                      << ", Target: " << min_cut_size << std::endl;
            break;
        }
    }
    
    std::cout << "Mincut enforcement completed. Total edges added: " << total_edges_added 
              << ", Degree corrected: " << degree_corrected_edges 
              << " (" << (total_edges_added > 0 ? (degree_corrected_edges * 100 / total_edges_added) : 0) 
              << "%)" << std::endl;
}

#endif // ENFORCE_MINCUT_WITH_BUDGET_H
