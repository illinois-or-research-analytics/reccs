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
    auto cluster_available_nodes = degree_aware_stages::get_cluster_available_nodes(task);
    
    std::cout << "Cluster has " << cluster_available_nodes.size() 
              << " nodes with available degree budget" << std::endl;

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
        std::set<std::pair<uint32_t, uint32_t>> edges_to_add;
        
        // Create sets for faster lookup
        std::set<unsigned int> light_set(light_partition.begin(), light_partition.end());
        std::set<unsigned int> heavy_set(heavy_partition.begin(), heavy_partition.end());
        
        // Use statics.h for efficient edge existence checking
        auto existing_edges = statics::compute_existing_edges(g);
        auto edge_exists = statics::create_edge_exists_checker(existing_edges);
        
        // Categorize nodes by partition and available budget
        struct PartitionNodes {
            std::vector<uint32_t> available_budget;
            std::vector<uint32_t> no_budget;
        };
        
        PartitionNodes light_nodes, heavy_nodes;
        
        // Categorize light partition nodes
        for (unsigned int node : light_partition) {
            uint64_t node_id = g.id_map[node];
            if (cluster_available_nodes.count(node_id) > 0 && 
                task.degree_manager->get_available_degree(node_id) > 0) {
                light_nodes.available_budget.push_back(node);
            } else {
                light_nodes.no_budget.push_back(node);
            }
        }
        
        // Categorize heavy partition nodes
        for (unsigned int node : heavy_partition) {
            uint64_t node_id = g.id_map[node];
            if (cluster_available_nodes.count(node_id) > 0 && 
                task.degree_manager->get_available_degree(node_id) > 0) {
                heavy_nodes.available_budget.push_back(node);
            } else {
                heavy_nodes.no_budget.push_back(node);
            }
        }
        
        std::cout << "Light partition: " << light_nodes.available_budget.size() 
                  << " with budget, " << light_nodes.no_budget.size() << " without" << std::endl;
        std::cout << "Heavy partition: " << heavy_nodes.available_budget.size() 
                  << " with budget, " << heavy_nodes.no_budget.size() << " without" << std::endl;
        
        // Generate candidate cross-partition edges with priority
        struct EdgeCandidate {
            uint32_t u, v;
            int priority; // Higher is better: 2=both have budget, 1=one has budget, 0=none have budget
            
            bool operator<(const EdgeCandidate& other) const {
                if (priority != other.priority) return priority > other.priority;
                if (u != other.u) return u < other.u;
                return v < other.v;
            }
        };
        
        std::vector<EdgeCandidate> candidate_edges;
        
        // Helper to add candidates between two node lists
        auto add_candidates = [&](const std::vector<uint32_t>& list1, 
                                 const std::vector<uint32_t>& list2, 
                                 int priority) {
            for (uint32_t u : list1) {
                for (uint32_t v : list2) {
                    if (!edge_exists(u, v)) {
                        uint32_t smaller = std::min(u, v);
                        uint32_t larger = std::max(u, v);
                        candidate_edges.push_back({smaller, larger, priority});
                    }
                }
            }
        };
        
        // Priority 2: Both nodes have available budget
        add_candidates(light_nodes.available_budget, heavy_nodes.available_budget, 2);
        
        // Priority 1: One node has available budget
        add_candidates(light_nodes.available_budget, heavy_nodes.no_budget, 1);
        add_candidates(light_nodes.no_budget, heavy_nodes.available_budget, 1);
        
        // Priority 0: Neither node has available budget
        add_candidates(light_nodes.no_budget, heavy_nodes.no_budget, 0);
        
        if (candidate_edges.empty()) {
            std::cout << "TERMINATION: No more edges can be added between partitions. "
                      << "All possible cross-partition edges already exist." << std::endl;
            std::cout << "Maximum achievable cut size with current partitioning: " << current_cut_size << std::endl;
            return;
        }
        
        // Sort by priority (highest first)
        std::sort(candidate_edges.begin(), candidate_edges.end());
        
        std::cout << "Found " << candidate_edges.size() << " possible new cross-partition edges:" << std::endl;
        std::cout << "  Priority 2 (both have budget): " 
                  << std::count_if(candidate_edges.begin(), candidate_edges.end(), 
                                  [](const EdgeCandidate& e) { return e.priority == 2; }) << std::endl;
        std::cout << "  Priority 1 (one has budget): " 
                  << std::count_if(candidate_edges.begin(), candidate_edges.end(), 
                                  [](const EdgeCandidate& e) { return e.priority == 1; }) << std::endl;
        std::cout << "  Priority 0 (none have budget): " 
                  << std::count_if(candidate_edges.begin(), candidate_edges.end(), 
                                  [](const EdgeCandidate& e) { return e.priority == 0; }) << std::endl;
        
        // Select edges to add (up to the number needed), preferring higher priority
        uint32_t edges_to_add_count = std::min(edges_needed, static_cast<uint32_t>(candidate_edges.size()));
        
        // Within each priority level, shuffle for randomness
        auto priority2_end = std::partition(candidate_edges.begin(), candidate_edges.end(),
                                           [](const EdgeCandidate& e) { return e.priority == 2; });
        auto priority1_end = std::partition(priority2_end, candidate_edges.end(),
                                           [](const EdgeCandidate& e) { return e.priority == 1; });
        
        std::shuffle(candidate_edges.begin(), priority2_end, gen);
        std::shuffle(priority2_end, priority1_end, gen);
        std::shuffle(priority1_end, candidate_edges.end(), gen);
        
        // Add the top priority edges
        for (uint32_t i = 0; i < edges_to_add_count; i++) {
            uint32_t u = candidate_edges[i].u;
            uint32_t v = candidate_edges[i].v;
            
            edges_to_add.insert({u, v});
            total_edges_added++;
            
            // Track if this edge uses available degree budget
            if (candidate_edges[i].priority > 0) {
                degree_corrected_edges++;
            }
        }
        
        // Convert to vector for batch addition
        std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
        
        if (edges_vector.empty()) {
            std::cout << "ERROR: No edges selected for addition." << std::endl;
            return;
        }
        
        std::cout << "Adding " << edges_vector.size() << " new cross-partition edges..." << std::endl;
        
        // Print breakdown of what we're adding
        size_t priority2_added = 0, priority1_added = 0, priority0_added = 0;
        for (uint32_t i = 0; i < edges_to_add_count; i++) {
            if (candidate_edges[i].priority == 2) priority2_added++;
            else if (candidate_edges[i].priority == 1) priority1_added++;
            else priority0_added++;
        }
        
        std::cout << "  Adding: " << priority2_added << " priority-2, " 
                  << priority1_added << " priority-1, " 
                  << priority0_added << " priority-0 edges" << std::endl;
        
        // Add the edges to the graph
        add_edges_batch(g, edges_vector);
        
        // Consume degree budgets for edges that were added
        for (const auto& edge : edges_vector) {
            uint64_t u_id = g.id_map[edge.first];
            uint64_t v_id = g.id_map[edge.second];
            
            // Try to consume budget (thread-safe)
            bool u_consumed = task.degree_manager->try_consume_degree(u_id, 1);
            bool v_consumed = task.degree_manager->try_consume_degree(v_id, 1);
            
            // Update cluster available nodes if needed
            if (u_consumed) {
                cluster_available_nodes.erase(u_id);
            }
            if (v_consumed) {
                cluster_available_nodes.erase(v_id);
            }
        }
        
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
