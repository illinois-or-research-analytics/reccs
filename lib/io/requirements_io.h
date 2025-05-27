#ifndef CONNECTIVITY_REQUIREMENTS_LOADER_H
#define CONNECTIVITY_REQUIREMENTS_LOADER_H

#include <unordered_map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <vector>
#include "mapped_file.h"

// Structure to hold cluster requirements
struct ClusterRequirements {
    std::string cluster_id;
    uint32_t num_nodes;
    uint32_t num_edges;
    uint32_t connectivity;
};

// Load connectivity requirements from CSV
class ConnectivityRequirementsLoader {
private:
    // Map from cluster ID to connectivity requirement
    std::unordered_map<std::string, uint32_t> cluster_to_connectivity;
    
    // Optional: Store full requirements if needed
    std::unordered_map<std::string, ClusterRequirements> cluster_requirements;
    
    // Parse a line from the CSV
    bool parse_line(const char* line_start, const char* line_end, int line_num) {
        // Skip empty lines
        if (line_start == line_end) return true;
        
        // Find commas
        const char* comma1 = std::strchr(line_start, ',');
        if (!comma1 || comma1 >= line_end) {
            std::cerr << "Error: Missing first comma on line " << line_num << std::endl;
            return false;
        }
        
        const char* comma2 = std::strchr(comma1 + 1, ',');
        if (!comma2 || comma2 >= line_end) {
            std::cerr << "Error: Missing second comma on line " << line_num << std::endl;
            return false;
        }
        
        const char* comma3 = std::strchr(comma2 + 1, ',');
        if (!comma3 || comma3 >= line_end) {
            std::cerr << "Error: Missing third comma on line " << line_num << std::endl;
            return false;
        }
        
        // Extract fields
        std::string cluster_id(line_start, comma1 - line_start);
        std::string n_str(comma1 + 1, comma2 - comma1 - 1);
        std::string m_str(comma2 + 1, comma3 - comma2 - 1);
        std::string conn_str(comma3 + 1, line_end - comma3 - 1);
        
        // Remove trailing whitespace/newlines
        while (!conn_str.empty() && (conn_str.back() == '\r' || conn_str.back() == '\n' || 
                                     conn_str.back() == ' ' || conn_str.back() == '\t')) {
            conn_str.pop_back();
        }
        
        // Convert to numbers
        try {
            uint32_t n = std::stoul(n_str);
            uint32_t m = std::stoul(m_str);
            uint32_t connectivity = std::stoul(conn_str);
            
            // Store the mapping
            cluster_to_connectivity[cluster_id] = connectivity;
            
            // Optionally store full requirements
            cluster_requirements[cluster_id] = {cluster_id, n, m, connectivity};
            
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error converting numbers on line " << line_num 
                      << ": " << e.what() << std::endl;
            return false;
        }
    }
    
public:
    // Load requirements from CSV file using memory mapping
    bool load_from_csv(const std::string& filename, bool verbose = false) {
        MappedFile file;
        if (!file.open(filename)) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }
        
        const char* data = file.data();
        const char* end = data + file.size();
        
        // Find first newline (end of header)
        const char* line_start = data;
        const char* line_end = std::strchr(line_start, '\n');
        
        if (!line_end || line_end >= end) {
            std::cerr << "Error: File appears to be empty or have no newline" << std::endl;
            return false;
        }
        
        // Verify header
        std::string header(line_start, line_end - line_start);
        if (!header.empty() && header.back() == '\r') {
            header.pop_back();
        }
        
        if (header != "cluster,n,m,connectivity") {
            std::cerr << "Warning: Unexpected header format: " << header << std::endl;
        }
        
        // Parse data lines
        int line_num = 1;
        line_start = line_end + 1;
        
        while (line_start < end) {
            line_num++;
            
            // Find end of line
            line_end = std::strchr(line_start, '\n');
            if (!line_end) line_end = end;
            
            // Parse line
            if (line_end > line_start) {
                parse_line(line_start, line_end, line_num);
            }
            
            // Move to next line
            line_start = line_end + 1;
        }
        
        if (verbose) {
            std::cout << "Loaded " << cluster_to_connectivity.size() 
                      << " cluster connectivity requirements" << std::endl;
        }
        
        return true;
    }
    
    // Get connectivity requirement for a cluster
    uint32_t get_connectivity_requirement(const std::string& cluster_id) const {
        auto it = cluster_to_connectivity.find(cluster_id);
        if (it != cluster_to_connectivity.end()) {
            return it->second;
        }
        // Return 0 if cluster not found
        return 0;
    }
    
    // Check if cluster has requirements
    bool has_requirements(const std::string& cluster_id) const {
        return cluster_to_connectivity.count(cluster_id) > 0;
    }
    
    // Get all cluster IDs
    std::vector<std::string> get_all_cluster_ids() const {
        std::vector<std::string> ids;
        ids.reserve(cluster_to_connectivity.size());
        for (const auto& [id, _] : cluster_to_connectivity) {
            ids.push_back(id);
        }
        return ids;
    }
    
    // Get full requirements for a cluster (if stored)
    const ClusterRequirements* get_full_requirements(const std::string& cluster_id) const {
        auto it = cluster_requirements.find(cluster_id);
        if (it != cluster_requirements.end()) {
            return &it->second;
        }
        return nullptr;
    }
    
    // Get the mapping directly
    const std::unordered_map<std::string, uint32_t>& get_connectivity_map() const {
        return cluster_to_connectivity;
    }
    
    // Clear all loaded data
    void clear() {
        cluster_to_connectivity.clear();
        cluster_requirements.clear();
    }
    
    // Get statistics
    void print_statistics() const {
        if (cluster_to_connectivity.empty()) {
            std::cout << "No connectivity requirements loaded" << std::endl;
            return;
        }
        
        uint32_t min_conn = UINT32_MAX;
        uint32_t max_conn = 0;
        uint64_t sum_conn = 0;
        
        for (const auto& [_, conn] : cluster_to_connectivity) {
            min_conn = std::min(min_conn, conn);
            max_conn = std::max(max_conn, conn);
            sum_conn += conn;
        }
        
        std::cout << "Connectivity Requirements Statistics:" << std::endl;
        std::cout << "  Number of clusters: " << cluster_to_connectivity.size() << std::endl;
        std::cout << "  Min connectivity: " << min_conn << std::endl;
        std::cout << "  Max connectivity: " << max_conn << std::endl;
        std::cout << "  Average connectivity: " 
                  << (double)sum_conn / cluster_to_connectivity.size() << std::endl;
    }
};

// Example usage
void example_usage() {
    ConnectivityRequirementsLoader loader;
    
    // Load from CSV
    if (!loader.load_from_csv("cluster_requirements.csv", true)) {
        std::cerr << "Failed to load requirements" << std::endl;
        return;
    }
    
    // Print statistics
    loader.print_statistics();
    
    // Get requirement for specific cluster
    std::string cluster_id = "100";
    if (loader.has_requirements(cluster_id)) {
        uint32_t connectivity = loader.get_connectivity_requirement(cluster_id);
        std::cout << "Cluster " << cluster_id << " requires connectivity: " 
                  << connectivity << std::endl;
    }
    
    // Iterate through all requirements
    for (const auto& [cluster_id, connectivity] : loader.get_connectivity_map()) {
        std::cout << "Cluster " << cluster_id << " -> " << connectivity << std::endl;
    }
}

#endif // CONNECTIVITY_REQUIREMENTS_LOADER_H
