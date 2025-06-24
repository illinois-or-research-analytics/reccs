#ifndef ENFORCE_MINCUT_H
#define ENFORCE_MINCUT_H

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

void enforce_mincut(Graph& g, uint32_t min_cut_size) {
    std::cout << "Starting mincut enforcement on cluster "
              << g.id << ". Target minimum cut size: " << min_cut_size << std::endl;

    int iteration = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    
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
        
        // Count existing cross-partition edges
        uint32_t existing_cross_count = 0;
        for (unsigned int u : light_partition) {
            for (unsigned int v : heavy_partition) {
                if (edge_exists(u, v)) {
                    existing_cross_count++;
                }
            }
        }
        
        std::cout << "Found " << existing_cross_count << " existing cross-partition edges." << std::endl;
        
        // Generate candidate cross-partition edges using efficient edge checking
        std::vector<std::pair<uint32_t, uint32_t>> candidate_edges;
        for (unsigned int u : light_partition) {
            for (unsigned int v : heavy_partition) {
                // Use O(1) edge existence check from statics.h
                if (!edge_exists(u, v)) {
                    uint32_t smaller = std::min(u, v);
                    uint32_t larger = std::max(u, v);
                    candidate_edges.push_back({smaller, larger});
                }
            }
        }
        
        if (candidate_edges.empty()) {
            std::cout << "TERMINATION: No more edges can be added between partitions. "
                      << "All possible cross-partition edges already exist." << std::endl;
            std::cout << "Maximum achievable cut size with current partitioning: " << current_cut_size << std::endl;
            return;
        }
        
        std::cout << "Found " << candidate_edges.size() << " possible new cross-partition edges." << std::endl;
        
        // Randomly shuffle candidates for variety
        std::shuffle(candidate_edges.begin(), candidate_edges.end(), gen);
        
        // Select edges to add (up to the number needed)
        uint32_t edges_to_add_count = std::min(edges_needed, static_cast<uint32_t>(candidate_edges.size()));
        
        for (uint32_t i = 0; i < edges_to_add_count; i++) {
            edges_to_add.insert(candidate_edges[i]);
        }
        
        // Convert to vector for batch addition
        std::vector<std::pair<uint32_t, uint32_t>> edges_vector(edges_to_add.begin(), edges_to_add.end());
        
        if (edges_vector.empty()) {
            std::cout << "ERROR: No edges selected for addition." << std::endl;
            return;
        }
        
        std::cout << "Adding " << edges_vector.size() << " new cross-partition edges..." << std::endl;
        
        // Add the edges to the graph
        add_edges_batch(g, edges_vector);
        
        std::cout << "Added " << edges_to_add.size() << " edges in iteration " << iteration << std::endl;
    }
}

#endif // ENFORCE_MINCUT_H
